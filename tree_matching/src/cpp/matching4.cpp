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


/*--------------------------------------------------------------*
 *
 *                      CIRAD / Modelisation
 *
 *  File : matching3.C
 *
 *  Class : Matching
 *
 *  Description: Classe implementant les algorithmes d'alignement
 *               entre arborescences, algorithme classique de Zhang,
 *               algorithme minimisant la distance et le nombre de 
 *               composantes connexes sur un (ou deux) arbre(s)
 *
 *  Author: Pascal Ferraro
 *
 *  Date: 19/07/99
 *
 *--------------------------------------------------------------*/
#include "matching4.h"


  // -------------
  // Constructeur
  // -------------
MatchingMinimizeComponents::MatchingMinimizeComponents(TreeGraph& input,TreeGraph& reference,NodeCost& nodeDistance)
{
  T1=&input;
  T2=&reference;
  ND=&nodeDistance;
  _distances.make(*T1,*T2,nodeDistance);
  _nbInputConnectedComponents.make(*T1,*T2);
  _nbTreeConnectedComponents.make(*T1,*T2);
  _nbInputConnectedComponentsT1.make(*T1,*T2);
  _nbReferenceConnectedComponents.make(*T1,*T2);
  _nbReferenceConnectedComponentsT2.make(*T1,*T2);
 
  _InputRootTable.make(*T1,*T2);
  _ReferenceRootTable.make(*T1,*T2);
  _InputRootTableT1.make(*T1,*T2);
  _ReferenceRootTableT2.make(*T1,*T2);

  _sumNbCaseVector = CaseVector(3,0);
  _nbCaseVector = CaseVector(7,0);
// _choices est un tableau de listes retenant les tentatives successives alignements durant l'algo.
// c'est donc un tableau de |T1| lignes et |T2| colonnes initialise a 0
  _choices.resize(T1->getNbVertex(),T2->getNbVertex());
  _choicesT1.resize(T1->getNbVertex(),T2->getNbVertex());
  _choicesT2.resize(T1->getNbVertex(),T2->getNbVertex());
// constante qui va permettre de calculer l'alignement restreint
  _restrMapp.link(I_MAX(T1->getDegree(),T2->getDegree()),*_distances.getDistanceTable(),*_nbTreeConnectedComponents.getConnectedTable());
}
  // -------------
  // Destructeur
  // -------------
MatchingMinimizeComponents::~MatchingMinimizeComponents()
{
}

  // ----------------------------------------------------------------------------------
  // Calcule la distance entre les deux arbres T1[input_vertex] et T2[reference_vertex]
  // ----------------------------------------------------------------------------------
DistanceType MatchingMinimizeComponents::distanceBetweenTree(int input_vertex,int reference_vertex)
{
  // On stocke dans ni et nj le nombre d'enfants de T1[i] et T2[j]
  int ni=T1->getNbChild(input_vertex);
  int nj=T2->getNbChild(reference_vertex);
  DistanceType cost1,cost2,cost3,dist1,dist2;
  DistanceType min,MIN=2*MAXDIST;
  int im=-1,jm=-1,MTC=0;
  int imT1=-1,jmT1=-1,MTCT1=0;
  int imT2=-1,jmT2=-1,MTCT2=0;
  int i;
  int nbcon,NBCON1=MAXINT,NBCON2=MAXINT,NBCON=MAXINT;
  int nbconref,NBCONREF1=MAXINT,NBCONREF2=MAXINT,NBCONREF=MAXINT;
  int nbconT1,NBCON1T1=MAXINT,NBCON2T1=MAXINT,NBCT1=MAXINT;
  int nbconrefT2,NBCONREF1T2=MAXINT,NBCONREF2T2=MAXINT,NBCREFT2=MAXINT;


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
  nbconT1=0;
  nbconrefT2=0;
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
      nbconT1 = getNBCT1(input_child,reference_vertex);
      nbconrefT2 = getNBCRefT2(input_child,reference_vertex);

      // On conserve la plus petite distance
      if (dist1<min)
	{ 
	  min=dist1; 
	  im=input_child; 
	  imT1=input_child; 
	  imT2=input_child; 
	  NBCON1=nbcon;
	  NBCONREF1=nbconref;
	  NBCON1T1=nbconT1;
	  NBCONREF1T2=nbconrefT2;
	}

      // Si la distance est egale, on conserve celui qui minimise le nbre
      // de composantes connexes de T1+T2
      if (dist1==min)
	{
	  if (nbcon+nbconref<NBCON1+NBCONREF1)
	    {
	      im=input_child; 
	      NBCON1=nbcon;
	      NBCONREF1=nbconref;
	    }
	  if (nbconT1<NBCON1T1)
	    {
	      imT1=input_child; 
	      NBCON1T1=nbconT1;
	    }
	  if (nbconrefT2<NBCONREF1T2)
	    {
	      imT2=input_child; 
	      NBCONREF1T2=nbconrefT2;
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
      NBCT1=NBCON1T1;
      MTCT1=1;
      MTCT2=1;
      NBCREFT2=NBCONREF1T2;
    }
  if (cost1==MIN)
    {
      if (NBCON1+NBCONREF1<NBCON+NBCONREF)
	{
	  MTC=1;
	  NBCON=NBCON1;
	  NBCONREF=NBCONREF1;
	}
      if (NBCON1T1<NBCT1)
	{
	  MTCT1=1;
	  NBCT1=NBCON1T1;
	}
      if (NBCONREF1T2<NBCREFT2)
	{
	  MTCT2=1;
	  NBCREFT2=NBCONREF1T2;
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
      nbconT1=getNBCT1(input_vertex,reference_child);
      nbconrefT2=getNBCRefT2(input_vertex,reference_child);
      
      if (dist2<min)
	{ 
	  min=dist2;
	  jm=reference_child;
	  NBCON2=nbcon;
	  NBCONREF2=nbconref;
	  jmT1=reference_child;
	  NBCON2T1=nbconT1;
	  jmT2=reference_child;
	  NBCONREF2T2=nbconrefT2;
	}
      if (dist2==min)
	{
	  if (nbcon+nbconref<NBCONREF2+NBCON2)
	    {
	      jm=reference_child;
	      NBCON2=nbcon;
	      NBCONREF2=nbconref;
	    }
	  if (nbconT1<NBCON2T1)
	    {
	      jmT1=reference_child;
	      NBCON2T1=nbconT1;
	    }
	  if (nbconrefT2<NBCONREF2T2)
	    {
	      jmT2=reference_child;
	      NBCONREF2T2=nbconrefT2;
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
      MTCT1=2;
      NBCT1=NBCON2T1;
      MTCT2=2;
      NBCREFT2=NBCONREF2T2;
    }
  if (cost2==MIN)
    {

      if (NBCON2+NBCONREF2<NBCONREF+NBCON)
	{
	  MTC=2;
	  NBCON=NBCON2;
	  NBCONREF=NBCONREF2;
	}
      if (NBCON2T1<=NBCT1)
	{
	  MTCT1=2;
	  NBCT1=NBCON2T1;
	}
      if (NBCONREF2T2<NBCREFT2)
	{
	  MTCT2=2;
	  NBCREFT2=NBCONREF2T2;
	}
    }


//----------------------------------------------------------------------------------
//Case3 : We evaluate the matching between the input_forest and the reference_forest
// On evalue la mise en correspondance des arbres des deux forets issues de T1 et T2
//----------------------------------------------------------------------------------
// Le cout est celui de l'alignement des deux forets
// plus celui de l'echange de T1(i) en T2(j)


  cost3=getDBF(input_vertex,reference_vertex);
  cost3=cost3+_distances.getCCost(input_vertex,reference_vertex);

  nbcon=_nbInputConnectedComponents.getNBCF(input_vertex,reference_vertex)
    -_InputRootTable.getNBCF(input_vertex,reference_vertex)
    +1;

  nbconref=_nbReferenceConnectedComponents.getNBCF(input_vertex,reference_vertex)
    -_ReferenceRootTable.getNBCF(input_vertex,reference_vertex)
    +1;

  nbconT1=_nbInputConnectedComponentsT1.getNBCF(input_vertex,reference_vertex)
    -_InputRootTableT1.getNBCF(input_vertex,reference_vertex)
    +1;

 
  nbconrefT2=_nbReferenceConnectedComponentsT2.getNBCF(input_vertex,reference_vertex)
    -_ReferenceRootTableT2.getNBCF(input_vertex,reference_vertex)
    +1;

  // On conserve le cout s'il est inferieur au precedent
  if (cost3<MIN) 
    { 
      MIN=cost3; 
      MTC=3; 
      NBCON=nbcon;
      NBCONREF=nbconref;
      MTCT1=3; 
      NBCT1=nbconT1;
      MTCT2=3; 
      NBCREFT2=nbconrefT2;
    }
  
  if (cost3==MIN)
    {

      if (nbcon+nbconref<=NBCONREF+NBCON)
	{
	  MTC=3;
	  NBCON=nbcon;
	  NBCONREF=nbconref;
	}
      if (nbconT1<=NBCT1)
	{
	  MTCT1=3;
	  NBCT1=nbconT1;
	}
      if (nbconrefT2<=NBCREFT2)
	{
	  MTCT2=3;
	  NBCREFT2=nbconrefT2;
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


  // Maintient des listes pour la minimisation dans les deux arbres 
  switch (MTCT1)
    {
    case 1 :{
      _choicesT1.putFirst(input_vertex,reference_vertex,imT1); 
      _choicesT1.putLast(input_vertex,reference_vertex,-1);
      _InputRootTableT1.putNBC(input_vertex,reference_vertex,0);
    };break;
    case 2 :{
      _choicesT1.putFirst(input_vertex,reference_vertex,jmT1);
      _choicesT1.putLast(input_vertex,reference_vertex,M(input_vertex,jmT1));
      _InputRootTableT1.putNBC(input_vertex,reference_vertex,1);
    };break;
    case 3 :{
      _choicesT1.putFirst(input_vertex,reference_vertex,-1);
      _choicesT1.putLast(input_vertex,reference_vertex,reference_vertex);
      _InputRootTableT1.putNBC(input_vertex,reference_vertex,1);
    };break;
    default : 	assert(0);break;
    }


  // Maintient des listes pour la minimisation dans l'arbre T2
  switch (MTCT2)
    {
    case 1 :{
      _choicesT2.putFirst(input_vertex,reference_vertex,imT2); 
      _choicesT2.putLast(input_vertex,reference_vertex,-1);
      _ReferenceRootTableT2.putNBC(input_vertex,reference_vertex,1);
    };break;
    case 2 :{
      _choicesT2.putFirst(input_vertex,reference_vertex,jmT2);
      _choicesT2.putLast(input_vertex,reference_vertex,M(input_vertex,jmT2));
      _ReferenceRootTableT2.putNBC(input_vertex,reference_vertex,0);
    };break;
    case 3 :{
      _choicesT2.putFirst(input_vertex,reference_vertex,-1);
      _choicesT2.putLast(input_vertex,reference_vertex,reference_vertex);
      _ReferenceRootTableT2.putNBC(input_vertex,reference_vertex,1);
    };break;
    default : 	assert(0);break;
    }
  _choices.putFirst(input_vertex,reference_vertex,MTC);
  _nbInputConnectedComponents.putNBC(input_vertex,reference_vertex,NBCON);
  _nbReferenceConnectedComponents.putNBC(input_vertex,reference_vertex,NBCONREF);
  _nbTreeConnectedComponents.putNBC(input_vertex,reference_vertex,NBCON+NBCONREF);
  

  _choicesT1.putFirst(input_vertex,reference_vertex,MTCT1);
  _nbInputConnectedComponentsT1.putNBC(input_vertex,reference_vertex,NBCT1);

  _choicesT2.putFirst(input_vertex,reference_vertex,MTCT2);
  _nbReferenceConnectedComponentsT2.putNBC(input_vertex,reference_vertex,NBCREFT2);

  _distances.putDBT(input_vertex,reference_vertex,MIN);

  /*if ((NBCT1!=NBCON)||(NBCREFT2!=NBCONREF))
    {
      cout<<" CONTRE EXEMPLE !!!!!!! "<<endl;
      cout<<"racines des deux arbres :"<<endl;
      (T1->getNode(input_vertex))->print();
      cout<<T1->getNbChild(input_vertex)<<endl;
      (T2->getNode(reference_vertex))->print();
      cout<<T2->getNbChild(reference_vertex)<<endl;
      cout<<"minimisation sur T1 :"<<NBCT1<<endl;
      cout<<"minimisation sur T2 :"<<NBCREFT2<<endl;
      cout<<"minimisation globale sur T1 et T2 : "<<NBCON<<" -  "<<NBCONREF<<endl;
      cout<<"cost1 ="<< cost1<<endl;
      cout<<"cost2 ="<< cost2<<endl;
      cout<<"cost3 ="<< cost3<<endl;
      cout<<MTC<<" - "<<MTCT1<<" - "<<MTCT2<<endl;


      }*/

  // On calcule les differents cas rencontres

  if ((cost1==cost2) && (cost1==cost3))
    { 
      _nbCaseVector[0]++;
    }
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
DistanceType MatchingMinimizeComponents::distanceBetweenForest(int input_vertex,int reference_vertex)
{
// ni et nj representent le nombre de forets a comparees
  int ni=T1->getNbChild(input_vertex);
  int nj=T2->getNbChild(reference_vertex);
  DistanceType cost1,cost2,cost3,dist1,dist2;
  DistanceType min,DIST;
  int im=-1,jm=-1,MFC=0;
  int i;
  int cf,CF1=MAXINT,CF2=MAXINT,CF3=MAXINT,CF=MAXINT; 
  int sa,SA2=0,SA=0;
  int sar,SAR1=0,SAR=0;
//cf est le nbre de composantes connexes de la foret initiale
  // il est egal a la somme des composantes connexes de
  //chaque arbre
  int cfref,CFREF1=MAXINT,CFREF2=MAXINT,CFREF3=MAXINT,CFREF=MAXINT;


  int imT1=-1,jmT1=-1,MFCT1=0;
  int cfT1,CF1T1=MAXINT,CF2T1=MAXINT,CF3T1=MAXINT,CFT1=MAXINT; 
  int saT1,SA2T1=0,SAT1=0;
  int imT2=-1,jmT2=-1,MFCT2=0;
  int sarT2,SAR1T2=0,SART2=0;
//cf est le nbre de composantes connexes de la foret initiale
  // il est egal a la somme des composantes connexes de
  //chaque arbre
  int cfrefT2,CFREF1T2=MAXINT,CFREF2T2=MAXINT,CFREF3T2=MAXINT,CFREFT2=MAXINT;
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

      cf=_nbInputConnectedComponents.getNBCF(input_child,reference_vertex);
      cfref=_nbReferenceConnectedComponents.getNBCF(input_child,reference_vertex);
      sar = _ReferenceRootTable.getNBCF(input_child,reference_vertex);

      cfT1=_nbInputConnectedComponentsT1.getNBCF(input_child,reference_vertex);
      cfrefT2=_nbReferenceConnectedComponentsT2.getNBCF(input_child,reference_vertex);
      sarT2 = _ReferenceRootTableT2.getNBCF(input_child,reference_vertex);



      if (dist1<min) 
	{ 
	  min=dist1;

	  im=input_child;
	  CF1=cf;
	  CFREF1=cfref;
	  SAR1=sar;

	  imT1=input_child;
	  CF1T1=cfT1;

	  imT2=input_child;
	  CFREF1T2=cfrefT2;
	  SAR1T2=sarT2;
	}
      if (dist1==min)
	{
	  if (cf+cfref-sar<CF1+CFREF1-SAR1)
	    {
	      im=input_child;
	      CF1=cf;
	      CFREF1=cfref;
	      SAR1=sar;
	    }
	  if (cfT1<CF1T1)
	    {
	      imT1=input_child;
	      CF1T1=cfT1;
	    }
	  if (cfrefT2-sarT2<CFREF1T2-SAR1T2)
	    {
	      imT2=input_child;
	      CFREF1T2=cfrefT2;
	      SAR1T2=sarT2;
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
      SA=0;
      _InputRootTable.putNBCF(input_vertex,reference_vertex,0);
      _ReferenceRootTable.putNBCF(input_vertex,reference_vertex,SAR1);


      MFCT1=1;
      CFT1=CF1T1;
      SAT1=0;
      _InputRootTableT1.putNBCF(input_vertex,reference_vertex,0);

      MFCT2=1;
      CFREFT2=CFREF1T2;
      SART2=SAR1T2;
      _ReferenceRootTableT2.putNBCF(input_vertex,reference_vertex,SAR1T2);
   }

  if (cost1==DIST) 
    {
      if (CF1+CFREF1-SAR1<CF+CFREF-SA-SAR)
	{
	  MFC=1;
	  CF=CF1;
	  CFREF=CFREF1;
	  SA=0;
	  SAR=SAR1;
	  _InputRootTable.putNBCF(input_vertex,reference_vertex,0);
	  _ReferenceRootTable.putNBCF(input_vertex,reference_vertex,SAR1);
	}

      if (CF1T1<CFT1-SAT1)
	{
	  MFCT1=1;
	  CFT1=CF1T1;
	  SAT1=0;
	  _InputRootTableT1.putNBCF(input_vertex,reference_vertex,0);
	}

      if (CFREF1T2-SAR1T2<CFREFT2-SART2)
	{
	  MFCT2=1;
	  CFREFT2=CFREF1T2;
	  SART2=SAR1T2;
	  _ReferenceRootTableT2.putNBCF(input_vertex,reference_vertex,SAR1T2);
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

      cf=_nbInputConnectedComponents.getNBCF(input_vertex,reference_child);
      cfref=_nbReferenceConnectedComponents.getNBCF(input_vertex,reference_child);
      sa=_InputRootTable.getNBCF(input_vertex,reference_child);


      cfT1=_nbInputConnectedComponentsT1.getNBCF(input_vertex,reference_child);
      saT1=_InputRootTableT1.getNBCF(input_vertex,reference_child);

      cfrefT2=_nbReferenceConnectedComponentsT2.getNBCF(input_vertex,reference_child);



      if (dist2<min) 
	{ 
	  min=dist2;
	  jm=reference_child;
	  CF2=cf;
	  CFREF2=cfref;
	  SA2=sa;

	  jmT1=reference_child;
	  CF2T1=cfT1;
	  SA2T1=saT1;

	  jmT2=reference_child;
	  CFREF2T2=cfrefT2;

	}
      if (dist2==min)
	{
	  if (cf+cfref-sa<CF2+CFREF2-SA2)
	    {
	      jm=reference_child;
	      CF2=cf;
	      CFREF2=cfref;
	      SA2=sa;
	    }

	  if (cfT1-saT1<CF2T1-SA2T1)
	    {
	      jmT1=reference_child;
	      CF2T1=cfT1;
	      SA2T1=saT1;
	    }

	  if (cfrefT2<CFREF2T2)
	    {
	      jmT2=reference_child;
	      CFREF2T2=cfrefT2;
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
      _InputRootTable.putNBCF(input_vertex,reference_vertex,SA2);
      _ReferenceRootTable.putNBCF(input_vertex,reference_vertex,0);
      SA=SA2;
      SAR=0;

      MFCT1=2;
      CFT1=CF2T1;
      _InputRootTableT1.putNBCF(input_vertex,reference_vertex,SA2T1);
      SAT1=SA2T1;

      MFCT2=2;
      CFREFT2=CFREF2T2;
      _ReferenceRootTableT2.putNBCF(input_vertex,reference_vertex,0);
      SART2=0;


    }
  if (cost2==DIST) 
    {
      if (CF2+CFREF2-SA2<CF+CFREF-SA-SAR)
	{
	  MFC=2;
	  CF=CF2;
	  CFREF=CFREF2;
	  SA=SA2;
	  SAR=0;
	  _InputRootTable.putNBCF(input_vertex,reference_vertex,SA2);
	  _ReferenceRootTable.putNBCF(input_vertex,reference_vertex,0);
	}


      if (CF2T1-SA2T1<CFT1-SAT1)
	{
	  MFCT1=2;
	  CFT1=CF2T1;
	  SAT1=SA2T1;
	  _InputRootTableT1.putNBCF(input_vertex,reference_vertex,SA2T1);
	}

      if (CFREF2T2<CFREFT2-SART2)
	{
	  MFCT2=2;
	  CFREFT2=CFREF2T2;
	  SART2=0;
	  _ReferenceRootTableT2.putNBCF(input_vertex,reference_vertex,0);
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


  int nbInputMappedChildT1=0;
  int nbReferenceMappedChildT2=0;



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

      CF3T1=0;
      CFREF3T2=0;
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

	  CF3T1=0;
	  CFREF3T2=0;

	}
      else
	{
	  //BOTH FOREST ARE NOT EMPTY_TREE
// A retricted mapping must be calculated
// Sinon on resout le probleme de flot maximum de cout minimum
	  cost3=_restrMapp.minCostFlow(_restrMappList);
	  //cost3=_restrMapp.minCostFlowWithComponents(_restrMappList);
	  cf=0;
	  cfref=0;
	  cfT1=0;
	  cfrefT2=0;
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




		  cfT1=cfT1+getNBCT1(input_child,_restrMapp.who(_restrMappList[i]));
		  nbInputMappedChildT1=nbInputMappedChildT1 
		    + _InputRootTableT1.getNBC(input_child,_restrMapp.who(_restrMappList[i]));

		  cfrefT2=cfrefT2+getNBCRefT2(input_child,_restrMapp.who(_restrMappList[i]));
		  nbReferenceMappedChildT2=nbReferenceMappedChildT2 
		    + _ReferenceRootTableT2.getNBC(input_child,_restrMapp.who(_restrMappList[i]));

		}
	    }
	  CF3=cf;
	  CFREF3=cfref;

	  CFREF3T2=cfrefT2;
	  CF3T1=cfT1;
	}
    }

  if (cost3<DIST)
    { 
      DIST=cost3; 
      MFC=3;
      CF=CF3;
      CFREF=CFREF3;
      _InputRootTable.putNBCF(input_vertex,reference_vertex,nbInputMappedChild);
      _ReferenceRootTable.putNBCF(input_vertex,reference_vertex,nbReferenceMappedChild);
      
      MFCT1=3;
      MFCT2=3;
      CFT1=CF3T1;
      CFREFT2=CFREF3T2;
      _InputRootTableT1.putNBCF(input_vertex,reference_vertex,nbInputMappedChildT1);
      _ReferenceRootTableT2.putNBCF(input_vertex,reference_vertex,nbReferenceMappedChildT2);

   }
  if (cost3==DIST)
    { 
      if (CF3+CFREF3-nbInputMappedChild-nbReferenceMappedChild<=CF+CFREF-SA-SAR)
	{
	  MFC=3;
	  CF=CF3;
	  CFREF=CFREF3;
	  _InputRootTable.putNBCF(input_vertex,reference_vertex,nbInputMappedChild);
	  _ReferenceRootTable.putNBCF(input_vertex,reference_vertex,nbReferenceMappedChild);
	}
      if (CF3T1-nbInputMappedChildT1<=CFT1-SAT1)
	{
	  MFCT1=3;
	  CFT1=CF3T1;
	  _InputRootTableT1.putNBCF(input_vertex,reference_vertex,nbInputMappedChildT1);
	}
      if (CFREF3T2-nbReferenceMappedChildT2<=CFREFT2-SART2)
	{ 
	  MFCT2=3;
	  CFREFT2=CFREF3T2;
	  _ReferenceRootTableT2.putNBCF(input_vertex,reference_vertex,nbReferenceMappedChildT2);
	}
      
    }

//-------------------------------
//We maintain the matching lists
// On maintient les listes d'alignement
//-------------------------------
  _choices.createList(input_vertex,reference_vertex);
  _choices.putFirst(input_vertex,reference_vertex,MFC);

  _nbInputConnectedComponents.putNBCF(input_vertex,reference_vertex,CF);
  _nbReferenceConnectedComponents.putNBCF(input_vertex,reference_vertex,CFREF);
	


  _choicesT1.createList(input_vertex,reference_vertex);
  _choicesT1.putFirst(input_vertex,reference_vertex,MFCT1);

  _nbInputConnectedComponentsT1.putNBCF(input_vertex,reference_vertex,CFT1);
	
  _choicesT2.createList(input_vertex,reference_vertex);
  _choicesT2.putFirst(input_vertex,reference_vertex,MFCT2);

  _nbReferenceConnectedComponentsT2.putNBCF(input_vertex,reference_vertex,CFREFT2);
	
  switch(MFC)
    {
    case 1 : _choices.putLast(input_vertex,reference_vertex,im);break;
    case 2 : _choices.putLast(input_vertex,reference_vertex,jm);break;
    case 3 : {
      for (int i=1;i<=T1->getNbChild(input_vertex);i++)
	{
	  _choices.putLast(input_vertex,reference_vertex,_restrMapp.who(_restrMappList[i]));
	}
    } break;
    default : 	break;
    }
	
  switch(MFCT1)
    {
    case 1 : _choicesT1.putLast(input_vertex,reference_vertex,imT1);break;
    case 2 : _choicesT1.putLast(input_vertex,reference_vertex,jmT1);break;
    case 3 : {
      for (int i=1;i<=T1->getNbChild(input_vertex);i++)
	{
	  _choicesT1.putLast(input_vertex,reference_vertex,_restrMapp.who(_restrMappList[i]));
	}
    } break;
    default : 	break;
    }
  switch(MFCT2)
    {
    case 1 : _choicesT2.putLast(input_vertex,reference_vertex,imT2);break;
    case 2 : _choicesT2.putLast(input_vertex,reference_vertex,jmT2);break;
    case 3 : {
      for (int i=1;i<=T1->getNbChild(input_vertex);i++)
	{
	  _choicesT2.putLast(input_vertex,reference_vertex,_restrMapp.who(_restrMappList[i]));
	}
    } break;
    default : 	break;
    }
    
  delete (NodeList*) input_list;
  delete (NodeList*) reference_list;
  _distances.putDBF(input_vertex,reference_vertex,DIST);
  return(DIST);

}



DistanceType MatchingMinimizeComponents::match()
{
  DistanceType D=0;
  if ((!T1->isNull())&&(!T2->isNull()))
    {
      for (int input_vertex=T1->getNbVertex()-1;input_vertex>=0;input_vertex--)
	{
	  _distances.openDistancesVector(input_vertex);
	  _InputRootTable.openIntVector(input_vertex);
	  _ReferenceRootTable.openIntVector(input_vertex);
	  _nbInputConnectedComponents.openIntVector(input_vertex);
	  _nbTreeConnectedComponents.openIntVector(input_vertex);
	  _nbReferenceConnectedComponents.openIntVector(input_vertex);

	  _InputRootTableT1.openIntVector(input_vertex);
	  _nbInputConnectedComponentsT1.openIntVector(input_vertex);

	  _ReferenceRootTableT2.openIntVector(input_vertex);
	  _nbReferenceConnectedComponentsT2.openIntVector(input_vertex);

	  //On initialise le tableau _InputRootTable a False, aucun noeud n'a d'images
	  //au debut de l'algo.
	  for (int j=0;j<T2->getNbVertex();j++)
	    {
	      _InputRootTable.putNBC(input_vertex,j,0);
	      _ReferenceRootTable.putNBC(input_vertex,j,0);
	      _InputRootTableT1.putNBC(input_vertex,j,0);
	      _ReferenceRootTableT2.putNBC(input_vertex,j,0);
	    }

	  for (int reference_vertex=T2->getNbVertex()-1;reference_vertex>=0;reference_vertex--)
	    {
	      distanceBetweenForest(input_vertex,reference_vertex);
	      distanceBetweenTree(input_vertex,reference_vertex);
	    }
	  
	  for (int i=1;i<=T1->getNbChild(input_vertex);i++)
	    {
	      _distances.closeDistancesVector(T1->child(input_vertex,i));
	      _nbInputConnectedComponents.closeIntVector(T1->child(input_vertex,i));
	      _nbTreeConnectedComponents.closeIntVector(T1->child(input_vertex,i));
	      _InputRootTable.closeIntVector(T1->child(input_vertex,i));
	      _nbReferenceConnectedComponents.closeIntVector(T1->child(input_vertex,i));
	      _ReferenceRootTable.closeIntVector(T1->child(input_vertex,i));
	      _nbInputConnectedComponentsT1.closeIntVector(T1->child(input_vertex,i));
	      _InputRootTableT1.closeIntVector(T1->child(input_vertex,i));

	      _nbReferenceConnectedComponentsT2.closeIntVector(T1->child(input_vertex,i));
	      _ReferenceRootTableT2.closeIntVector(T1->child(input_vertex,i));
	    }
	}
      
      D=getDBT(0,0);
      cout<<"minimisation sur T1 :"<<getNBCT1(0,0)<<std::endl;
      cout<<"minimisation sur T2 :"<<getNBCRefT2(0,0)<<std::endl;
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
	cout<<endl;*/
      cout<<"nbre de composantes connexes sur T1 : "<<getNBC(0,0)<<endl;
      cout<<"nbre de composantes connexes sur T2 : "<<getNBCRef(0,0)<<endl;

      
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

int MatchingMinimizeComponents::getNBC(int input_vertex,int reference_vertex) const
{
  return(_nbInputConnectedComponents.getNBC(input_vertex,reference_vertex));
}

int MatchingMinimizeComponents::getNBCRef(int input_vertex,int reference_vertex) const
{
  return(_nbReferenceConnectedComponents.getNBC(input_vertex,reference_vertex));
}



int MatchingMinimizeComponents::getNBCT1(int input_vertex,int reference_vertex) const
{
  return(_nbInputConnectedComponentsT1.getNBC(input_vertex,reference_vertex));
}

int MatchingMinimizeComponents::getNBCRefT2(int input_vertex,int reference_vertex) const
{
  return(_nbReferenceConnectedComponentsT2.getNBC(input_vertex,reference_vertex));
}


