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


#include "matching_with_complex.h"
#include "dec_matrix.h"


  // -------------
  // Constructeur
  // -------------
MatchingWithComplex::MatchingWithComplex(TreeGraph& input,TreeGraph& reference,NodeCost& nodeDistance)
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
  _restrMapp_v_w.link(I_MAX(T1->getDegree(),T2->getDegree()),*_restrDistances_v_w.getDistanceTable());
  _restrMapp_l_l.link(I_MAX(T1->getDegree(),T2->getDegree()),*_restrDistances_l_l.getDistanceTable());
}



// -------------
// Destructeur
  // -------------
MatchingWithComplex::~MatchingWithComplex()
{
}

// ----------------------------------------------------------------------------------
// Calcule la distance entre les deux arbres T1[input_vertex] et T2[reference_vertex]
// tel que le complexe de v a une image mais pas celui de w
// ----------------------------------------------------------------------------------
DistanceType MatchingWithComplex::distanceBetweenTree(int input_vertex,int reference_vertex)
{
  // On stocke dans ni et nj le nombre d'enfants de T1[i] et T2[j]
  int ni=T1->getNbChild(input_vertex);
  int nj=T2->getNbChild(reference_vertex);
  
  DistanceType cost,dist1,dist2,dist;
  DistanceType MIN_v_l=MAXDIST,MIN_l_w=MAXDIST;
  DistanceType MIN_l_l=MAXDIST,MIN_v_w=MAXDIST,MIN=MAXDIST;
  DistanceType min,min_v_l,min_l_l,min_l_w,min_v_w;
  int im=0,jm=0,jm_v_l=0,im_l_w=0,MTC=0;
  int im_v_w=0,jm_v_w=0;
  int MTC_l_w=0,MTC_v_l=0,MTC_l_l=0,MTC_v_w=0;
  int mtc=0,mtc_l_w=0,mtc_v_l=0,mtc_l_l=0,mtc_v_w=0;
  int i;

  // Evaluation de _d_l_w
	
  cost=getDBT(input_vertex,EMPTY_TREE);
  min     = MAXDIST; // pour le calcul de _d
  min_l_w = MAXDIST; // pour le calcul de _d_l_w
  min_v_w = MAXDIST; // pour le calcul de _d_v_w
  
  for (i=1;i<=ni;i++)
    {
      int input_child=T1->child(input_vertex,i);
      
      dist  = _distances.getDBT(input_child,reference_vertex)-getDBT(input_child,EMPTY_TREE); 
      dist1 = _d_l_w.getDBT(input_child,reference_vertex)-getDBT(input_child,EMPTY_TREE);
      dist2 = _d_v_w.getDBT(input_child,reference_vertex)-getDBT(input_child,EMPTY_TREE);

      if (dist<min)
	{
	  min = dist;
	  im  = input_child;
	  mtc = 1;
	}
      if (dist1<min_l_w)
	{
	  min_l_w = dist1;
	  im_l_w  = input_child;
	  mtc_l_w = 1;
	}
      if (T1->childIsInComplex(input_vertex,i))
	{
	  if (dist2<min_v_w)
	    {
	      min_v_w = dist2;
	      im_v_w  = input_child;
	      mtc_v_w = 1;
	    }
	}
      else
	{
	  if (dist2<min_l_w)
	    {
	      min_l_w = dist2;
	      im_l_w  = input_child;
	      mtc_l_w = 2;
	    }  
	}
    }
  min     = min + cost;
  min_l_w = min_l_w + cost;
  min_v_w = min_v_w+cost;

  if (min_l_w<MIN_l_w)
    {
      MIN_l_w = min_l_w;
      MTC_l_w = mtc_l_w;                  // c'est l'ancien cas 1
    }
  if (min_v_w<MIN_v_w)
    {
      MIN_v_w = min_v_w;
      MTC_v_w = mtc_v_w;                  // c'est l'ancien cas 1
    }
  if (min<MIN)
    {
      MIN = min;
      MTC = mtc;                  // c'est l'ancien cas 1
    }
  

  cost=getDBT(EMPTY_TREE,reference_vertex);
  min     = MAXDIST; // pour le calcul de _d
  min_v_l = MAXDIST; // pour le calcul de _d_v_l
  min_v_w = MAXDIST; // pour le calcul de _d_v_w
  for (i=1;i<=nj;i++)
    {
      int reference_child=T2->child(reference_vertex,i);
      
      dist1 = _d_v_l.getDBT(input_vertex,reference_child)-getDBT(EMPTY_TREE,reference_child);
      dist2 = _d_v_w.getDBT(input_vertex,reference_child)-getDBT(EMPTY_TREE,reference_child);
      dist  = _distances.getDBT(input_vertex,reference_child)-getDBT(EMPTY_TREE,reference_child);

      if (dist<min)
	{
	  min = dist;
	  jm  = reference_child;
	  mtc = 2;
	}
      
      if (dist1<min_v_l)
	{
	  min_v_l = dist1;
	  jm_v_l  = reference_child;
	  mtc_v_l = 1;
	}
      
  
      if (T2->childIsInComplex(reference_vertex,i))
	{
	  if (dist2<min_v_w)
	    {
	      min_v_w = dist2;
	      jm_v_w  = reference_child;
	      mtc_v_w = 2;
	    }
	}
      else
	{
	  if (dist2<min_v_l)
	    {
	      min_v_l = dist2;
	      jm_v_l  = reference_child;
	      mtc_v_l = 2;
	    }
	}
    }
  min     = min + cost;
  min_v_l = min_v_l+cost;
  min_v_w = min_v_w+cost;

  if (min_v_l<MIN_v_l)
    {
      MIN_v_l = min_v_l;
      MTC_v_l = mtc_v_l;                  // c'est l'ancien cas 2
    }
  if (min_v_w<MIN_v_w)
    {
      MIN_v_w = min_v_w;
      MTC_v_w = mtc_v_w;                  // c'est l'ancien cas 2
    }
  if (min <MIN )
    {
      MIN  = min ;
      MTC  = mtc;                  // c'est l'ancien cas 2
    }
  
  min_v_l = _d_v_l.getDBF(input_vertex,reference_vertex)+_distances.getICost(reference_vertex)
    +_distances.getDCost(input_vertex);
  mtc_v_l = 3;

  min_l_w = _d_l_w.getDBF(input_vertex,reference_vertex)+_distances.getDCost(input_vertex)
    +_distances.getICost(reference_vertex);
  mtc_l_w = 3;

  
  min_l_l = _d_l_l.getDBF(input_vertex,reference_vertex)+_distances.getDCost(input_vertex)
    +_distances.getICost(reference_vertex);
  mtc_l_w = 1;

  if (_d_v_w.getDBF(input_vertex,reference_vertex)<=_d_l_l.getDBF(input_vertex,reference_vertex))
    {
      min_v_w = _d_v_w.getDBF(input_vertex,reference_vertex);
      mtc_v_w = 3;
    }
  
  else
    {
      min_v_w = _d_l_l.getDBF(input_vertex,reference_vertex);
      mtc_v_w = 4;
    }
  

  min_v_w = min_v_w + _distances.getCCost(input_vertex,reference_vertex);

  if (min_v_l<MIN_v_l)
    {
      MIN_v_l = min_v_l;
      MTC_v_l = mtc_v_l;
    }
  if (min_l_w<MIN_l_w)
    {
      MIN_l_w = min_l_w;
      MTC_l_w = mtc_l_w;
    }
  if (min_v_w<MIN_v_w)
    {
      MIN_v_w = min_v_w;
      MTC_v_w = mtc_v_w;
    }
  if (min_l_l<MIN_l_l)
    {
      MIN_l_l = min_l_l;
      MTC_l_l = mtc_l_l;
    }
  if (min_l_w <MIN )
    {
      MIN  = min_l_w ;
      MTC  = 5;
    }
  if (min_v_l <MIN )
    {
      MIN  = min_v_l ;
      MTC  = 6;
    }
  if (min_v_w <MIN )
    {
      MIN  = min_v_w ;
      MTC  = mtc_v_w;
    }

      
      
  switch (MTC)
    {
    case 1 :{
      _choices.putFirst(input_vertex,reference_vertex,im); 
      _choices.putLast(input_vertex,reference_vertex,-1); 
    }
    break;
    case 2 :{
      _choices.putFirst(input_vertex,reference_vertex,jm);
      _choices.putLast(input_vertex,reference_vertex,M(input_vertex,jm));
    }
    break;
    case 3 : case 4 : {
      _choices.putFirst(input_vertex,reference_vertex,-1);
      _choices.putLast(input_vertex,reference_vertex,reference_vertex);
    }
    break;
    case 5 : case 6 :{
      _choices.putFirst(input_vertex,reference_vertex,-1);
      _choices.putLast(input_vertex,reference_vertex,-1);
    }
    break;
    default: 	assert(0);break;
    }
  _choices.putFirst(input_vertex,reference_vertex,MTC);
  // On range dans le tableau des distances, la distance entre les arbres de racines input_vertex et
  // reference_vertex.
  _distances.putDBT(input_vertex,reference_vertex,MIN);
  
  switch (MTC_v_l)
    {
    case 1 : case 2 :{
      _choices_v_l.putFirst(input_vertex,reference_vertex,jm_v_l);
      _choices_v_l.putLast(input_vertex,reference_vertex,M_v_l(input_vertex,jm_v_l));
    }
    break;
    case 3 :{
      _choices_v_l.putFirst(input_vertex,reference_vertex,-1);
      _choices_v_l.putLast(input_vertex,reference_vertex,-1);
    }
    break;
    default: 	assert(0);break;
    }
  _choices_v_l.putFirst(input_vertex,reference_vertex,MTC_v_l);
  _d_v_l.putDBT(input_vertex,reference_vertex,MIN_v_l);
  
  switch (MTC_l_w)
    {
    case 1 : case 2 :{
      _choices_l_w.putFirst(input_vertex,reference_vertex,im_l_w); 
      _choices_l_w.putLast(input_vertex,reference_vertex,-1); 
    }
    break;
    case 3 :{
      _choices_l_w.putFirst(input_vertex,reference_vertex,-1);
      _choices_l_w.putLast(input_vertex,reference_vertex,-1);
    }
    break;
    default: 	assert(0);break;
    }
  _choices_l_w.putFirst(input_vertex,reference_vertex,MTC_l_w);
  _d_l_w.putDBT(input_vertex,reference_vertex,MIN_l_w);

  
  switch (MTC_v_w)
    {
    case 1 :{
      _choices_v_w.putFirst(input_vertex,reference_vertex,im_v_w); 
      _choices_v_w.putLast(input_vertex,reference_vertex,-1); 
    }
    break;
    case 2 :{
      _choices_v_w.putFirst(input_vertex,reference_vertex,jm_v_w);
      _choices_v_w.putLast(input_vertex,reference_vertex,M_v_w(input_vertex,jm_v_w));
    }
   
    break;

    case 3 : case 4 :{
      _choices_v_w.putFirst(input_vertex,reference_vertex,-1);
      _choices_v_w.putLast(input_vertex,reference_vertex,reference_vertex);
    }
    break;
    default: 	assert(0);break;
    }
  _choices_v_w.putFirst(input_vertex,reference_vertex,MTC_v_w);
  _d_v_w.putDBT(input_vertex,reference_vertex,MIN_v_w);

  _choices_l_l.putFirst(input_vertex,reference_vertex,-1);
  _choices_l_l.putLast(input_vertex,reference_vertex,-1);
  _choices_l_l.putFirst(input_vertex,reference_vertex,MTC_l_l);
  _d_l_l.putDBT(input_vertex,reference_vertex,MIN_l_l);
    
   return(MIN);
	
}


// -------------------------------
// Calcule la distance entre deux
// forets
// -------------------------------
DistanceType MatchingWithComplex::distanceBetweenForest(int input_vertex,int reference_vertex)
{
  // On stocke dans ni et nj le nombre d'enfants de T1[i] et T2[j]
  int ni=T1->getNbChild(input_vertex);
  int nj=T2->getNbChild(reference_vertex);
  
  DistanceType indel_forests,cost1,cost2,cost3,cost,dist1,dist2,dist3,dist4;
  DistanceType MIN=MAXDIST,MIN_v_l=MAXDIST,MIN_l_w=MAXDIST;
  DistanceType MIN_l_l=MAXDIST,MIN_v_w=MAXDIST;
  DistanceType min;
  DistanceType min_v_l,min_l_w,min_l_l,min_v_w;
  int im=0,jm=0,im_v_l=0,jm_v_l=0,im_l_w=0,jm_l_w=0;
  int im_v_w=0,jm_v_w=0,im_l_l=0,jm_l_l=0,imv=0,jmv=0,imw=0,jmw=0;
  int MFC_l_w=0,MFC_v_l=0,MFC_l_l=0,MFC_v_w=0,MFC = 0;
  int mtc_l_w=0,mtc_v_l=0,mtc_l_l=0,mtc_v_w=0;
  int i,j;

	
  cost=getDBF(input_vertex,EMPTY_TREE);
  min_l_w = MAXDIST; 
  min_v_l = MAXDIST; 
  min_v_w = MAXDIST; 
  min_l_l = MAXDIST; 
  
  for (i=1;i<=ni;i++)
    {
      int input_child=T1->child(input_vertex,i);

      dist1 = _d_l_w.getDBF(input_child,reference_vertex)-getDBF(input_child,EMPTY_TREE);
      dist2 = _d_v_l.getDBF(input_child,reference_vertex)-getDBF(input_child,EMPTY_TREE);
      dist3 = _d_v_w.getDBF(input_child,reference_vertex)-getDBF(input_child,EMPTY_TREE);
      dist4 = _d_l_l.getDBF(input_child,reference_vertex)-getDBF(input_child,EMPTY_TREE);
      
      if (dist1<min_l_w)
	{
	  min_l_w = dist1;
	  im_l_w = input_child;
	  mtc_l_w = 1;
	}
      if (dist4<min_l_l)
	{
	  min_l_l = dist4;
	  im_l_l = input_child;
	  mtc_l_l = 1;
	}

     if (T1->childIsInComplex(input_vertex,i))
	{
	  if (dist2<min_v_l)
	    {
	      min_v_l = dist2;
	      im_v_l = input_child;
	      mtc_v_l = 1;
	    }
	  if (dist3<min_v_w)                  
	    {
	      min_v_w = dist3;
	      im_v_w = input_child;
	      mtc_v_w = 1;
	    }
	}
     else
	{
	  if (dist3<min_l_w)
	    {
	      min_l_w = dist3;
	      im_l_w  = input_child;
	      mtc_l_w = 2;
	    }
	  if (dist2<min_l_l)
	    {
	      min_l_l = dist2;
	      im_l_l  = input_child;
	      mtc_l_l = 2;
	    }
	}
    }


  min_l_w = min_l_w+cost;
  min_v_l = min_v_l+cost;
  min_v_w = min_v_w+cost;
  min_l_l = min_l_l+cost;
    


  if (min_l_w<MIN_l_w)
    {
      MIN_l_w = min_l_w;
      MFC_l_w = mtc_l_w;                 
    }
    

  if (min_v_l<MIN_v_l)
    {
      MIN_v_l = min_v_l;
      MFC_v_l = mtc_v_l;                 
    }
  if (min_v_w<MIN_v_w)
    {
      MIN_v_w = min_v_w;
      MFC_v_w = mtc_v_w;                 
    }
  if (min_l_l<MIN_l_l)
    {
      MIN_l_l = min_l_l;
      MFC_l_l = mtc_l_l;                 
    }
 

  cost=getDBT(EMPTY_TREE,reference_vertex);
  min_v_l = MAXDIST;
  min_l_w = MAXDIST; 
  min_v_w = MAXDIST; 
  min_l_l = MAXDIST; 
  
  for (i=1;i<=nj;i++)
    {
      int reference_child=T2->child(reference_vertex,i);
      
      dist1=_d_v_l.getDBF(input_vertex,reference_child)-getDBF(EMPTY_TREE,reference_child);
      dist2=_d_l_w.getDBF(input_vertex,reference_child)-getDBF(EMPTY_TREE,reference_child);
      dist3=_d_v_w.getDBF(input_vertex,reference_child)-getDBF(EMPTY_TREE,reference_child);
      dist4=_d_l_l.getDBF(input_vertex,reference_child)-getDBF(EMPTY_TREE,reference_child);


      if (dist1<min_v_l)
	{
	  min_v_l = dist1;
	  jm_v_l  = reference_child;
	  mtc_v_l = 2;
	}

       if (dist4<min_l_l)
	{
	  min_l_l = dist4;
	  jm_l_l  = reference_child;
	  mtc_l_l = 3;
	}
     
      if (T2->childIsInComplex(reference_vertex,i))
	{
	  if (dist2<min_l_w)
	    {
	      min_l_w = dist2;
	      jm_l_w  = reference_child;
	      mtc_l_w = 3;
	    }
	  if (dist3<min_v_w)
	    {
	      min_v_w = dist3;
	      jm_v_w  = reference_child;
	      mtc_v_w = 2;
	    }
	}
      else
	{
	  if (dist3<min_v_l)
	    {
	      min_v_l = dist3;
	      jm_v_l  = reference_child;
	      mtc_v_l = 3;
	    }
	  if (dist2<min_l_l)
	    {
	      min_l_l = dist2;
	      mtc_l_l = 4;
	      jm_l_l  = reference_child;
	    }
	}
    }


  min_l_w = min_l_w+cost;
  min_v_l = min_v_l+cost;
  min_v_w = min_v_w+cost;
  min_l_l = min_l_l+cost;
    


  if (min_l_w<MIN_l_w)
    {
      MIN_l_w = min_l_w;
      MFC_l_w = mtc_l_w;                 
    }
    

  if (min_v_l<MIN_v_l)
    {
      MIN_v_l = min_v_l;
      MFC_v_l = mtc_v_l;                 
    }
  if (min_v_w<MIN_v_w)
    {
      MIN_v_w = min_v_w;
      MFC_v_w = mtc_v_w;                 
    }
  if (min_l_l<MIN_l_l)
    {
      MIN_l_l = min_l_l;
      MFC_l_l = mtc_l_l;                 
    }

 
  
  // On fabrique le graphe de flot necessaire a la resolution du probleme

  // Calcul de l'alignement restreint dans le cas ou les complexes des racines n'ont pas d'image :
  // _d_l_l

  NodeList* input_list=new NodeList();
  NodeList* reference_list=new NodeList();
  

  IntVector choice(ni+1);
  NEW_MAT(int,choice_l_l,ni+1,nj+1)
  NEW_MAT(int,choice_v_w,ni+1,nj+1)
  // int choice_l_l[ni+1][nj+1];
  // int choice_v_w[ni+1][nj+1];



  for (int s1=1;s1<=ni;s1++) 
    input_list->push_back(T1->child(input_vertex,s1));

  for (int s2=1;s2<=nj;s2++) 
    reference_list->push_back(T2->child(reference_vertex,s2)); 

  _restrMapp_l_l.make(*input_list,*reference_list);
  _restrMappList_l_l.resize(ni+nj+3,EMPTY_NODE);


  _restrMapp_v_w.make(*input_list,*reference_list);
  _restrMappList_v_w.resize(ni+nj+3,EMPTY_NODE);

    
  if (ni==0) 
    { 
      _restrMappList_l_l[1]=2;
      _restrMappList_v_w[1]=2;
      for (i=1;i<=nj;i++)
	{
	  _restrMappList_l_l[i+1]=1;
	  _restrMappList_v_w[i+1]=1;
	}	  
      min_l_l=getDBF(EMPTY_TREE,reference_vertex);
      min_v_l=getDBF(EMPTY_TREE,reference_vertex);
      min_l_w=getDBF(EMPTY_TREE,reference_vertex);
      min_v_w=getDBF(EMPTY_TREE,reference_vertex);
    }
  else
    {
      if (nj==0) 
	{ 
	  _restrMappList_l_l[2]=1;
	  _restrMappList_v_w[2]=1;
	  for (i=1;i<=ni;i++)
	    {
	    _restrMappList_l_l[i]=ni+1;
	    _restrMappList_v_w[i]=ni+1;	
	    }
	  min_l_l=getDBF(input_vertex,EMPTY_TREE); 
	  min_v_l=getDBF(input_vertex,EMPTY_TREE); 
	  min_l_w=getDBF(input_vertex,EMPTY_TREE); 
	  min_v_w=getDBF(input_vertex,EMPTY_TREE); 
	}
      else
	{
	  // mise a jour de la matrice de distance dans le cas du calcull de _d_l_l
	  dist1 = MAXDIST;
	  dist2 = MAXDIST;
	  indel_forests = getDBF(input_vertex,EMPTY_TREE) + getDBF(EMPTY_TREE,reference_vertex);
		int s1;
	  for (s1=1;s1<=ni;s1++) 
	    {
	      int input_child = T1->child(input_vertex,s1);
	      for (int s2=1;s2<=nj;s2++) 
		{
		  int reference_child = T2->child(reference_vertex,s2);
		  min = _distances.getDBT(input_child,reference_child);
		  min_l_w = _d_l_w.getDBT(input_child,reference_child);
		  min_v_l = _d_v_l.getDBT(input_child,reference_child);
		  min_l_l = _d_l_l.getDBT(input_child,reference_child);
		  min_v_w = _d_v_w.getDBT(input_child,reference_child);
		  cost1 = min_v_l+indel_forests-getDBT(input_child,EMPTY_TREE)-getDBF(EMPTY_TREE,reference_child);
		  cost2 = min_l_w+indel_forests-getDBT(input_child,EMPTY_TREE)-getDBF(EMPTY_TREE,reference_child);
		  cost3 = min_v_w+indel_forests-getDBT(input_child,EMPTY_TREE)-getDBF(EMPTY_TREE,reference_child);
		  if (T1->childIsInComplex(input_vertex,s1))
		    {
		      if (T2->childIsInComplex(reference_vertex,s2))
			{
			  // Mise a jour du graphe de flot dans le cas ou pi(v) et pi(w) n'ont pas d'images
			  _restrDistances_l_l.putDBT(input_child,reference_child,min_l_l);
			  choice_l_l[s1][s2]= 4 ;
			  // Mise a jour du graphe de flot dans le cas ou pi(v) et pi(w) ontune image
			  _restrDistances_v_w.putDBT(input_child,reference_child,min_v_w);
			  choice_v_w[s1][s2]= 3 ;
			  // Mise a jour du graphe de flot dans le cas où pi(v) a une image mais pas  pi(w)
			  if (cost1<dist1)
			    {
			      dist1 = cost1;
			      imv = input_child;
			      jmv = reference_child;
			      mtc_v_l = 4;
			    }			      
			  // Mise a jour du graphe de flot dans le cas où pi(w) a une image mais pas  pi(v)
			  if (cost2<dist2)
			    {
			      dist2 = cost2;
			      imw = input_child;
			      jmw = reference_child;
			      mtc_l_w = 4;
			    }
			}
		      else
			{
			  if (min_l_w<min_l_l)
			    {
			      _restrDistances_l_l.putDBT(input_child,reference_child,min_l_w);
			      _restrDistances_v_w.putDBT(input_child,reference_child,min_l_w);
			      choice_l_l[s1][s2]= 2 ;
			      choice_v_w[s1][s2]= 2 ;
			    }
			  else
			    {
			      _restrDistances_l_l.putDBT(input_child,reference_child,min_l_l);
			      _restrDistances_v_w.putDBT(input_child,reference_child,min_l_l);
			      choice_l_l[s1][s2]= 4 ;
			      choice_v_w[s1][s2]= 4 ;
			    }	
			  if (cost1<dist1)
			    {
			      dist1 = cost1;
			      imv = input_child;
			      jmv = reference_child;
			      mtc_v_l = 4;
			    }			      
			  if (cost3<dist1)
			    {
			      dist1 = cost3;
			      imv = input_child;
			      jmv = reference_child;
			      mtc_v_l = 5;
			    }
			}
		    }
		  else		    
 		    {
		      if (T2->childIsInComplex(reference_vertex,s2))
			{
			  if (min_v_l<min_l_l)
			    {
			      _restrDistances_l_l.putDBT(input_child,reference_child,min_v_l);
			      _restrDistances_v_w.putDBT(input_child,reference_child,min_v_l);
			      choice_l_l[s1][s2]= 1 ;
			      choice_v_w[s1][s2]= 1 ;
			    }
			  else
			    {
			      _restrDistances_l_l.putDBT(input_child,reference_child,min_l_l);
			      _restrDistances_v_w.putDBT(input_child,reference_child,min_l_l);
			      choice_l_l[s1][s2]= 4 ;
			      choice_v_w[s1][s2]= 4 ;
			    }	
			  if (cost2<dist2)
			    {
			      dist2 = cost2;
			      imw = input_child;
			      jmw = reference_child;
			      mtc_l_w = 4;
			    }
			  if (cost3<dist2)
			    {
			      dist2 = cost3;
			      imw = input_child;
			      jmw = reference_child;
			      mtc_l_w = 5;
			    }
			}
		      else
			{
			  _restrDistances_l_l.putDBT(input_child,reference_child,min);
			  _restrDistances_v_w.putDBT(input_child,reference_child,min);
			  choice_l_l[s1][s2]= 5 ;
			  choice_v_w[s1][s2]= 5 ;
			}
		    }
		}
	    }	
	  min_v_l = dist1;
	  min_l_w = dist2;
	  min_l_l=_restrMapp_l_l.minCostFlow(_restrMappList_l_l);
	  min_v_w = MAXDIST; 
	  for ( s1=1;s1<=ni;s1++) 
	    {
	      _restrDistances = _restrDistances_v_w;
	      int input_child = T1->child(input_vertex,s1);
	      if (T1->childIsInComplex(input_vertex,s1))
		{
		  for (int s2=1;s2<=nj;s2++) 
		    {
		      int reference_child = T2->child(reference_vertex,s2);
		      if (!T2->childIsInComplex(reference_vertex,s2))
			{
			  _restrDistances.putDBT(input_child,reference_child,MAXDIST);
			  
			}
		    }
		  
		  min = _restrMapp_v_w.minCostFlow(_restrMappList_v_w);
		  if (min<min_v_w)
		    {
		      min_v_w = min;
		      for (int i=1;i<=T1->getNbChild(input_vertex);i++)
			{
			  choice[i] = _restrMapp_v_w.who(_restrMappList_v_w[i]);
			}
		    }
		}
	    }
	}
    }
  


  if (min_l_w<MIN_l_w)
    {
      MIN_l_w = min_l_w;
      MFC_l_w = mtc_l_w;                 
    }
    

  if (min_v_l<MIN_v_l)
    {
      MIN_v_l = min_v_l;
      MFC_v_l = mtc_v_l;                 
    }
  if (min_v_w<MIN_v_w)
    {
      MIN_v_w = min_v_w;
      MFC_v_w = 3;                 
    }
  if (min_l_l<MIN_l_l)
    {
      MIN_l_l = min_l_l;
      MFC_l_l = 5;                 
    }

  if (MIN_l_w<MIN)
    {
      MIN = MIN_l_w;
      im = im_l_w;
      MFC = MFC_l_w;
    }
  if (MIN_v_l<MIN)
    {
      MIN = MIN_v_l;
      im = im_v_l;
      MFC = MFC_v_l;
    }
  if (MIN_v_w<MIN)
    {
      MIN = MIN_v_w;
      im = im_v_w;
      MFC = MFC_v_w;
    }
  if (MIN_l_l<MIN)
    {
      MIN = MIN_l_l;
      im = im_l_l;
      MFC = MFC_l_l;
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
	    if (T1->child(input_vertex,i) == im)
	      _choices.putLast(input_vertex,reference_vertex,jm);
	    else
	      _choices.putLast(input_vertex,reference_vertex,-1);
	      
	  }
      }
      break;
    default : 	break;
    }

  _choices_v_l.createList(input_vertex,reference_vertex);
  _choices_v_l.putFirst(input_vertex,reference_vertex,MFC_v_l);

  switch(MFC_v_l)
    {
    case 1 : 
      {
	_choices_v_l.putLast(input_vertex,reference_vertex,im_v_l);
      }
      break;
    case 2 : case 3 :
      {
	_choices_v_l.putLast(input_vertex,reference_vertex,jm_v_l);
      }
      break;
    case 4 : case 5 : 
      {
	for (int i=1;i<=T1->getNbChild(input_vertex);i++)
	  {
	    if (T1->child(input_vertex,i) == imv)
	      _choices_v_l.putLast(input_vertex,reference_vertex,jmv);
	    else
	      _choices_v_l.putLast(input_vertex,reference_vertex,-1);
	      
	  }
      }
      break;
    default : 	break;
    }
  _choices_l_w.createList(input_vertex,reference_vertex);
  _choices_l_w.putFirst(input_vertex,reference_vertex,MFC_l_w);

  switch(MFC_l_w)
    {
    case 1 : case 2 : 
      {
	_choices_l_w.putLast(input_vertex,reference_vertex,im_l_w);
      }
      break;
    case 3 :
      {
	_choices_l_w.putLast(input_vertex,reference_vertex,jm_l_w);
      }
      break;
    case 4 : case 5 : 
	{
	  for (int i=1;i<=T1->getNbChild(input_vertex);i++)
	    {
	      if (T1->child(input_vertex,i) == imw)
		_choices_l_w.putLast(input_vertex,reference_vertex,jmw);
	      else
		_choices_l_w.putLast(input_vertex,reference_vertex,-1);
	    }
	}
	break;
    default : 	break;
    }

  _choices_l_l.createList(input_vertex,reference_vertex);
  _choices_l_l.putFirst(input_vertex,reference_vertex,MFC_l_l);
  switch(MFC_l_l)
    {
    case 1 : case 2 : 
      {
	_choices_l_l.putLast(input_vertex,reference_vertex,im_l_l);
      }
      break;
    case 3 : case 4 :
      {
	_choices_l_l.putLast(input_vertex,reference_vertex,jm_l_l);
      }
      break;
    case 5 : 
      {
	for (i=1;i<=T1->getNbChild(input_vertex);i++)
	  {
	    int reference_child = _restrMapp_l_l.who(_restrMappList_l_l[i]);
	    if (reference_child==-1)
	      {
		_choices_l_l.putLast(input_vertex,reference_vertex,-1);
	      }
	    else
	      {
		for (j=1;j<=T2->getNbChild(reference_vertex);j++)
		  {
		    if (T2->child(reference_vertex,j)==reference_child)
		      {
			_choices_l_l.putLast(input_vertex,reference_vertex,choice_l_l[i][j]);
		      }
		  }
	      }
	    _choices_l_l.putLast(input_vertex,reference_vertex,reference_child);
	  }
      }
      break;
    default : 	break;
    }


  _choices_v_w.createList(input_vertex,reference_vertex);
  _choices_v_w.putFirst(input_vertex,reference_vertex,MFC_v_w);
  
  switch(MFC_v_w)
    {
    case 1 : 
      {
	_choices_v_w.putLast(input_vertex,reference_vertex,im_v_w);
      }
      break;
    case 2 :
      {
	_choices_v_w.putLast(input_vertex,reference_vertex,jm_v_w);
      }
      break;
    case 3 : 
      {
	for (int i=1;i<=T1->getNbChild(input_vertex);i++)
	  {
	    if (choice[i]==-1)
	      {
		_choices_v_w.putLast(input_vertex,reference_vertex,-1);
	      }
	    else
	      {
		for (j=1;j<=T2->getNbChild(reference_vertex);j++)
		  {
		    if (T2->child(reference_vertex,j)==choice[i])
		      {
			_choices_v_w.putLast(input_vertex,reference_vertex,choice_v_w[i][j]);
		      }
		  }
	      }
	    _choices_v_w.putLast(input_vertex,reference_vertex,choice[i]);
	  }
      }
	break;
    default : 	break;
    }
  
  
  delete (NodeList*) input_list;
  delete (NodeList*) reference_list;
  DEL_MAT(choice_l_l,ni+1)
  DEL_MAT(choice_v_w,ni+1)
  _d_v_l.putDBF(input_vertex,reference_vertex,MIN_v_l);
  _d_l_w.putDBF(input_vertex,reference_vertex,MIN_l_w);
  _d_v_w.putDBF(input_vertex,reference_vertex,MIN_v_w);
  _d_l_l.putDBF(input_vertex,reference_vertex,MIN_l_l);
  // Les matrices de distance ne sont utilises que pour le calcul de l'alignement restreint
  // il  n'y a en fait pas de distance entre foret correspondante, il est neanmoins
  // necessaire de les initialisees
  _restrDistances_l_l.putDBF(input_vertex,reference_vertex,0.);
  _restrDistances_v_w.putDBF(input_vertex,reference_vertex,0.);
  return(MIN);
  
}


void MatchingWithComplex::getList(int input_vertex, int reference_vertex, Sequence* sequence)
{
  TreeList(input_vertex,reference_vertex,*sequence);
}



int MatchingWithComplex::Lat(ChoiceList* L, int vertex)
{
  ChoiceList::iterator begin;
  begin = L->begin();
  for (int i=0;i<vertex;i++)
    begin++;
  return(*begin);
}



void MatchingWithComplex::TreeList(int input_vertex,int reference_vertex,Sequence& sequence)
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
	case 3:
	  {
	    sequence.append(input_vertex,reference_vertex,_distances.getCCost(input_vertex,reference_vertex));
	    ForestList_v_w(input_vertex,reference_vertex,sequence);
	  }
	  break;
	case 4:
	  {
	    sequence.append(input_vertex,reference_vertex,_distances.getCCost(input_vertex,reference_vertex));
	    ForestList_l_l(input_vertex,reference_vertex,sequence);
	  }
	  break;
	case 5:
	  {
	    ForestList_v_l(input_vertex,reference_vertex,sequence);
	  }
	case 6:
	  {
	    ForestList_l_w(input_vertex,reference_vertex,sequence);
	  }
	default : break;
	}
    }
}


void MatchingWithComplex::TreeList_v_l(int input_vertex,int reference_vertex,Sequence& sequence)
{
  //cout<<"TreeList_v_l"<<input_vertex<<"  "<<reference_vertex<<endl;
  if ((!T1->isNull())&&(!T2->isNull()))
    {
      ChoiceList* L=_choices_v_l.getList(input_vertex,reference_vertex);
      int tree_choice=L->front();
      
      switch(tree_choice)
	{
	case 1: 
	  {
	    TreeList_v_l(Lat(L,1),reference_vertex,sequence);
	  }
	  break;
	case 2: 
	  {
	    TreeList_v_w(input_vertex,Lat(L,1),sequence);
	  }
	  break;
	case 3:
	  {
	    ForestList_v_l(input_vertex,reference_vertex,sequence);
	  }
	  break;
	default : break;
	}
    }
}

void MatchingWithComplex::TreeList_l_w(int input_vertex,int reference_vertex,Sequence& sequence)
{
  //cout<<"TreeList_l_w"<<input_vertex<<"  "<<reference_vertex<<endl;
  if ((!T1->isNull())&&(!T2->isNull()))
    {
      ChoiceList* L=_choices_l_w.getList(input_vertex,reference_vertex);
      int tree_choice=L->front();
      
      switch(tree_choice)
	{
	case 1: 
	  {
	    TreeList_l_w(Lat(L,1),reference_vertex,sequence);
	  }
	  break;
	case 2: 
	  {
	    TreeList_l_l(input_vertex,Lat(L,1),sequence);
	  }
	  break;
	case 3:
	  {
	    ForestList_l_w(input_vertex,reference_vertex,sequence);
	  }
	  break;
	default : break;
	}
    }
}
void MatchingWithComplex::TreeList_v_w(int input_vertex,int reference_vertex,Sequence& sequence)
{
  if ((!T1->isNull())&&(!T2->isNull()))
    {
      ChoiceList* L=_choices_v_w.getList(input_vertex,reference_vertex);
      int tree_choice=L->front();
      
      switch(tree_choice)
	{
	case 1: 
	  {
	    TreeList_v_w(Lat(L,1),reference_vertex,sequence);
	  }
	  break;
	case 2: 
	  {
	    TreeList_v_w(input_vertex,Lat(L,1),sequence);
	  }
	  break;
	case 3:
	  {
	    sequence.append(input_vertex,reference_vertex,_distances.getCCost(input_vertex,reference_vertex));
	    ForestList_v_w(input_vertex,reference_vertex,sequence);
	  }
	  break;
	case 4:
	  {
	    sequence.append(input_vertex,reference_vertex,_distances.getCCost(input_vertex,reference_vertex));
	    ForestList_l_l(input_vertex,reference_vertex,sequence);
	  }
	  break;
	default : break;
	}
    }
}
void MatchingWithComplex::TreeList_l_l(int input_vertex,int reference_vertex,Sequence& sequence)
{
  if ((!T1->isNull())&&(!T2->isNull()))
    {
      ChoiceList* L=_choices.getList(input_vertex,reference_vertex);
      int tree_choice=L->front();
      
      ForestList_l_l(input_vertex,reference_vertex,sequence);

    }
}

// Renvoie les matching listes entre forêts

void MatchingWithComplex::ForestList_v_l(int input_vertex,int reference_vertex,Sequence& sequence)
{
  //cout<<"ForestList_v_l"<<input_vertex<<"  "<<reference_vertex<<endl;
  ChoiceList* L=_choices_v_l.getList(input_vertex,reference_vertex);
  int forest_choice=Lat(L,2);
  
  switch(forest_choice)
    {
    case 1: 
      {
	ForestList_v_l(Lat(L,3),reference_vertex,sequence);
      }
      break;
    case 2: 
      {
	ForestList_v_l(input_vertex,Lat(L,3),sequence);
      }
      break;
    case 3: 
      {
	ForestList_v_w(input_vertex,Lat(L,3),sequence);
      }
      break;
    case 4: 
      {
	for (int i=1;i<=T1->getNbChild(input_vertex);i++)
	  {
	    int i_node=T1->child(input_vertex,i);
	    int r_node=Lat(L,2+i);
	    if (r_node!=-1) TreeList_v_l(i_node,r_node,sequence);
	  }
      }
      break;
    case 5: 
      {
	for (int i=1;i<=T1->getNbChild(input_vertex);i++)
	  {
	    int i_node=T1->child(input_vertex,i);
	    int r_node=Lat(L,2+i);
	    if (r_node!=-1) TreeList_v_w(i_node,r_node,sequence);
	  }
      }
      break;
    default : break;
    }
}


void MatchingWithComplex::ForestList_l_w(int input_vertex,int reference_vertex,Sequence& sequence)
{
  //cout<<"ForestList_l_w"<<input_vertex<<"  "<<reference_vertex<<endl;
  ChoiceList* L=_choices_l_w.getList(input_vertex,reference_vertex);
  int forest_choice=Lat(L,2);
  
  switch(forest_choice)
    {
    case 1: 
      {
	ForestList_l_w(Lat(L,3),reference_vertex,sequence);
      }
      break;
    case 2: 
      {
	ForestList_v_w(Lat(L,3),reference_vertex,sequence);
      }
      break;
    case 3: 
      {
	ForestList_l_w(input_vertex,Lat(L,3),sequence);
      }
      break;
    case 4: 
      {
	for (int i=1;i<=T1->getNbChild(input_vertex);i++)
	  {
	    int i_node=T1->child(input_vertex,i);
	    int r_node=Lat(L,2+i);
	    if (r_node!=-1) TreeList_l_w(i_node,r_node,sequence);
	  }
      }
      break;
    case 5: 
      {
	for (int i=1;i<=T1->getNbChild(input_vertex);i++)
	  {
	    int i_node=T1->child(input_vertex,i);
	    int r_node=Lat(L,2+i);
	    if (r_node!=-1) TreeList_v_w(i_node,r_node,sequence);
	  }
      }
      break;
    default : break;
    }
}
void MatchingWithComplex::ForestList_l_l(int input_vertex,int reference_vertex,Sequence& sequence)
{
  ChoiceList* L=_choices_l_l.getList(input_vertex,reference_vertex);
  int forest_choice=Lat(L,2);
  
  switch(forest_choice)
    {
    case 1: 
      {
	ForestList_l_l(Lat(L,3),reference_vertex,sequence);
      }
      break;
    case 2: 
      {
	ForestList_v_l(Lat(L,3),reference_vertex,sequence);
      }
      break;
    case 3: 
      {
	ForestList_l_l(input_vertex,Lat(L,3),sequence);
      }
      break;
    case 4: 
      {
	ForestList_l_w(input_vertex,Lat(L,3),sequence);
      }
      break;
    case 5: 
      {
	for (int i=1;i<=T1->getNbChild(input_vertex);i++)
	  {
	    int choice=Lat(L,2+2*i-1);
	    int i_node=T1->child(input_vertex,i);
	    int r_node=Lat(L,2+2*i);
	    switch(choice)
	      {
	      case 1: 
		{
		  if (r_node!=-1) TreeList_v_l(i_node,r_node,sequence);
		}
		break;
	      case 2: 
		{
		  if (r_node!=-1) TreeList_l_w(i_node,r_node,sequence);
		}
		break;
	      case 3: 
		{
		  if (r_node!=-1) TreeList_v_w(i_node,r_node,sequence);
		}
		break;
	      case 4: 
		{
		  if (r_node!=-1) TreeList_l_l(i_node,r_node,sequence);
		}
		break;
	      case 5: 
		{
		  if (r_node!=-1) TreeList(i_node,r_node,sequence);
		}
		break;
	      }
	  }
      }
      break;
    default : break;
    }
}


void MatchingWithComplex::ForestList_v_w(int input_vertex,int reference_vertex,Sequence& sequence)
{
  ChoiceList* L=_choices_v_w.getList(input_vertex,reference_vertex);
  int forest_choice=Lat(L,2);
  
  switch(forest_choice)
    {
    case 1: 
      {
	ForestList_v_l(Lat(L,3),reference_vertex,sequence);
      }
      break;
    case 2: 
      {
	ForestList_l_w(input_vertex,Lat(L,3),sequence);
      }
      break;
    case 3: 
      {
	for (int i=1;i<=T1->getNbChild(input_vertex);i++)
	  {
	    int choice=Lat(L,2+2*i-1);
	    int i_node=T1->child(input_vertex,i);
	    int r_node=Lat(L,2+2*i);
	    switch(choice)
	      {
	      case 1: 
		{
		  if (r_node!=-1) TreeList_v_l(i_node,r_node,sequence);
		}
		break;
	      case 2: 
		{
		  if (r_node!=-1) TreeList_l_w(i_node,r_node,sequence);
		}
		break;
	      case 3: 
		{
		  if (r_node!=-1) TreeList_v_w(i_node,r_node,sequence);
		}
		break;
	      case 4: 
		{
		  if (r_node!=-1) TreeList_l_l(i_node,r_node,sequence);
		}
		break;
	      case 5: 
		{
		  if (r_node!=-1) TreeList(i_node,r_node,sequence);
		}
		break;
	      }
	  }
      }
      break;
    default : break;
    }
}

DistanceType MatchingWithComplex::match()
{
  
  DistanceType D=0;
  if ((!T1->isNull())&&(!T2->isNull()))
    {
      for (int input_vertex=T1->getNbVertex()-1;input_vertex>=0;input_vertex--)
	{
	  _distances.openDistancesVector(input_vertex);
	  _d_v_l.openDistancesVector(input_vertex);
	  _d_l_w.openDistancesVector(input_vertex);
	  _d_v_w.openDistancesVector(input_vertex);
	  _d_l_l.openDistancesVector(input_vertex);
	  _restrDistances_v_w.openDistancesVector(input_vertex);
	  _restrDistances_l_l.openDistancesVector(input_vertex);
	  _restrDistances.openDistancesVector(input_vertex);
	  
	  for (int reference_vertex=T2->getNbVertex()-1;reference_vertex>=0;reference_vertex--)
	    {
	      distanceBetweenForest(input_vertex,reference_vertex);
	      distanceBetweenTree(input_vertex,reference_vertex);
	    }
	  
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

DistanceType MatchingWithComplex::getDBF(int input_vertex,int reference_vertex) const
{
  return(_distances.getDBF(input_vertex,reference_vertex));
}

DistanceType MatchingWithComplex::getDBT(int input_vertex,int reference_vertex) const
{
  return(_distances.getDBT(input_vertex,reference_vertex));
}


// renvoie le dernier element de la liste de la case node du tableau maintenant les listes d'alignement
int MatchingWithComplex::M(int i_node,int r_node)
{
  return(_choices.getList(i_node,r_node)->back());
}		

int MatchingWithComplex::M_v_l(int i_node,int r_node)
{
  return(_choices_v_l.getList(i_node,r_node)->back());
}		

int MatchingWithComplex::M_v_w(int i_node,int r_node)
{
  return(_choices_v_w.getList(i_node,r_node)->back());
}		

int MatchingWithComplex::M_l_w(int i_node,int r_node)
{
  return(_choices_l_w.getList(i_node,r_node)->back());
}		

int MatchingWithComplex::M_l_l(int i_node,int r_node)
{
  return(_choices_l_l.getList(i_node,r_node)->back());
}		





