#include "multiscale_matching.h"


  // -------------
  // Constructeur
  // -------------
MultiscaleMatching::MultiscaleMatching(TreeGraph& input,TreeGraph& reference,NodeCost& nodeDistance)
{
  T1=&input;
  T2=&reference;
  ND=&nodeDistance;
  _d.make(*T1,*T2,nodeDistance);
  _d_l_w.make(*T1,*T2,nodeDistance);
  _d_v_l.make(*T1,*T2,nodeDistance);
  _d_l_l_v_l.make(*T1,*T2,nodeDistance);
  _d_l_l_l_w.make(*T1,*T2,nodeDistance);
  _d_l_l_v_w.make(*T1,*T2,nodeDistance);
  _d_l_l_l_l.make(*T1,*T2,nodeDistance);
  _d_v_w_v_l.make(*T1,*T2,nodeDistance);
  _d_v_w_l_w.make(*T1,*T2,nodeDistance);
  _d_v_w_v_w.make(*T1,*T2,nodeDistance);
  _d_v_w_.make(*T1,*T2,nodeDistance);
  _d_v_w.make(*T1,*T2,nodeDistance);
  _d_l_l.make(*T1,*T2,nodeDistance);
  _restrDistances_v_w_.make(*T1,*T2,nodeDistance);
  _restrDistances.make(*T1,*T2,nodeDistance);
  _restrDistances_l_l_v_w.make(*T1,*T2,nodeDistance);
  _restrDistances_v_w_v_w.make(*T1,*T2,nodeDistance);
 // _choices est un tableau de listes retenant les tentatives successives alignements durant l'algo.
  // c'est donc un tableau de |T1| lignes et |T2| colonnes initialise a 0
  _choices.resize(T1->getNbVertex(),T2->getNbVertex());
  _choices_v_l.resize(T1->getNbVertex(),T2->getNbVertex());
  _choices_l_w.resize(T1->getNbVertex(),T2->getNbVertex());
  _choices_v_w_v_l.resize(T1->getNbVertex(),T2->getNbVertex());
  _choices_v_w_l_w.resize(T1->getNbVertex(),T2->getNbVertex());
  _choices_v_w_v_w.resize(T1->getNbVertex(),T2->getNbVertex());
  _choices_v_w_.resize(T1->getNbVertex(),T2->getNbVertex());
  _choices_v_w.resize(T1->getNbVertex(),T2->getNbVertex());
  _choices_l_l.resize(T1->getNbVertex(),T2->getNbVertex());
  _choices_l_l_v_l.resize(T1->getNbVertex(),T2->getNbVertex());
  _choices_l_l_l_w.resize(T1->getNbVertex(),T2->getNbVertex());
  _choices_l_l_v_w.resize(T1->getNbVertex(),T2->getNbVertex());
  _choices_l_l_l_l.resize(T1->getNbVertex(),T2->getNbVertex());
  // constante qui va permettre de calculer l'alignement restreint
  _restrMapp_v_w_.link(I_MAX(T1->getDegree(),T2->getDegree()),*_restrDistances_v_w_.getDistanceTable());
  _restrMapp.link(I_MAX(T1->getDegree(),T2->getDegree()),*_restrDistances.getDistanceTable());
  _restrMapp_l_l_v_w.link(I_MAX(T1->getDegree(),T2->getDegree()),*_restrDistances_l_l_v_w.getDistanceTable());
  _restrMapp_v_w_v_w.link(I_MAX(T1->getDegree(),T2->getDegree()),*_restrDistances_v_w_v_w.getDistanceTable());
}
// -------------
// Destructeur
  // -------------
MultiscaleMatching::~MultiscaleMatching()
{
}

// ----------------------------------------------------------------------------------
// Calcule la distance entre les deux arbres T1[input_vertex] et T2[reference_vertex]
// quand pi(v) n'a pas d'image et pi(w) a une image
// ----------------------------------------------------------------------------------
DistanceType MultiscaleMatching::distanceBetweenTree_v_l(int input_vertex,int reference_vertex)
{
  // On stocke dans ni et nj le nombre d'enfants de T1[i] et T2[j]
  int ni=T1->getNbChild(input_vertex);
  int nj=T2->getNbChild(reference_vertex);
  
  DistanceType cost2,cost3;
  DistanceType min,MIN=2*MAXDIST;
  int jm=0,MTC=0,mtc=0;
  int i;
  DistanceType dist_v_l,dist_v_w;

  min=MAXDIST;

  cost2=getDBT(EMPTY_TREE,reference_vertex);   
  for (i=1;i<=nj;i++)     
    {       
      int reference_child=T2->child(reference_vertex,i);
      dist_v_l=_d_v_l.getDBT(input_vertex,reference_child)-getDBT(EMPTY_TREE,reference_child);    
      dist_v_w=_d_v_w.getDBT(input_vertex,reference_child)-getDBT(EMPTY_TREE,reference_child);    
       if (dist_v_l<min) 
	{
	  min=dist_v_l;
	  jm=reference_child;
	  mtc = 1;
	}
       if (!(T2->childIsInComplex(reference_vertex,i)))
	{
	  if (dist_v_w<min)
	    {
	      min = dist_v_w;
	      jm  = reference_child;
	      mtc = 2;
	    }
	}
    }               
  cost2=cost2+min; 

  if (cost2<MIN) 
    {
      MIN=cost2; 
      MTC=mtc; 
    }


  cost3 = _d_v_l.getDBF(input_vertex,reference_vertex)+_d.getICost(reference_vertex)
    +_d.getDCost(input_vertex);
  if (cost3<MIN) 
    { 
      MIN=cost3;
      MTC=3; 
    }
  
  //-----------------------------------
  // We maintain the matching lists
  // mise a jour des listes d'alignement
  //-----------------------------------
  switch (MTC)
    {
    case 1 : case 2 :{
      _choices_v_l.putFirst(input_vertex,reference_vertex,jm);
      _choices_v_l.putLast(input_vertex,reference_vertex,M_v_l(input_vertex,jm));
    }
    break;
    case 3 :{
      _choices_v_l.putFirst(input_vertex,reference_vertex,-1);
      _choices_v_l.putLast(input_vertex,reference_vertex,-1);
    }
    break;
    default: 	assert(0);break;
    }
  _choices_v_l.putFirst(input_vertex,reference_vertex,MTC);
  _d_v_l.putDBT(input_vertex,reference_vertex,MIN);
  return(MIN);
}

// ----------------------------------------------------------------------------------
// Calcule la distance entre les deux arbres T1[input_vertex] et T2[reference_vertex]
// quand pi(w) n'a pas d'image et pi(v) a une image
// ----------------------------------------------------------------------------------
DistanceType MultiscaleMatching::distanceBetweenTree_l_w(int input_vertex,int reference_vertex)
{
  // On stocke dans ni et nj le nombre d'enfants de T1[i] et T2[j]
  int ni=T1->getNbChild(input_vertex);
  int nj=T2->getNbChild(reference_vertex);
  
  DistanceType cost1,cost3;
  DistanceType min,MIN=2*MAXDIST;
  int im=0,MTC=0, mtc =0;
  int i;
  DistanceType dist,dist_l_w,dist_v_w;

  min=MAXDIST;

  cost1=getDBT(input_vertex,EMPTY_TREE);   
  for (i=1;i<=ni;i++)
    {
      int input_child=T1->child(input_vertex,i);
      dist_l_w = _d_l_w.getDBT(input_child,reference_vertex)-getDBT(input_child,EMPTY_TREE);
      dist_v_w = _d_v_w.getDBT(input_child,reference_vertex)-getDBT(input_child,EMPTY_TREE);

      if (dist_l_w<min) 
	{ 
	  min=dist_l_w; 
	  im=input_child; 
	  mtc = 1;
	}
      if (!(T1->childIsInComplex(input_vertex,i)))
	{
	  if (dist_v_w<min)
	    {
	      min = dist_v_w;
	      im  = input_child;
	      mtc = 2;
	    }  
	}
    }
  cost1=cost1+min;
  // On conserve le cout minimum
  if (cost1<MIN)
    {
      MIN=cost1; 
      MTC=mtc; 
    }

  cost3 = _d_l_w.getDBF(input_vertex,reference_vertex)+_d.getICost(reference_vertex)
    +_d.getDCost(input_vertex);
  if (cost3<MIN) 
    { 
      MIN=cost3;
      MTC=3; 
    }
  
  //-----------------------------------
  // We maintain the matching lists
  // mise a jour des listes d'alignement
  //-----------------------------------
  switch (MTC)
    {
    case 1 : case 2 :{
      _choices_l_w.putFirst(input_vertex,reference_vertex,im); 
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
  _choices_l_w.putFirst(input_vertex,reference_vertex,MTC);
  _d_l_w.putDBT(input_vertex,reference_vertex,MIN);
  return(MIN);
}


// ----------------------------------------------------------------------------------
// Calcule la distance entre les deux arbres T1[input_vertex] et T2[reference_vertex]
// quand pi(v) n'a pas d'image et pi(w) a une image
// ----------------------------------------------------------------------------------
DistanceType MultiscaleMatching::distanceBetweenTree_v_w_v_l(int input_vertex,int reference_vertex)
{
  // On stocke dans ni et nj le nombre d'enfants de T1[i] et T2[j]
  int ni=T1->getNbChild(input_vertex);
  int nj=T2->getNbChild(reference_vertex);
  
  DistanceType cost1,cost2,cost3;
  DistanceType min,MIN=2*MAXDIST;
  int im=0,jm=0,MTC=0,mtc=0;
  int i;
  DistanceType dist,dist_v_w_v_l,dist_v_w;

  min=MAXDIST;

  cost1=getDBT(input_vertex,EMPTY_TREE);   
  for (i=1;i<=ni;i++)
    {
      int input_child=T1->child(input_vertex,i);
      dist_v_w_v_l = _d_v_w_v_l.getDBT(input_child,reference_vertex)-getDBT(input_child,EMPTY_TREE);

      if (T1->childIsInComplex(input_vertex,i))
	{
	  if (dist_v_w_v_l<min)
	    {
	      min = dist_v_w_v_l;
	      im  = input_child;
	      mtc = 1;
	    }  
	}
    }
  cost1=cost1+min;
  // On conserve le cout minimum
  if (cost1<MIN)
    {
      MIN=cost1; 
      MTC=mtc; 
    }

  min=MAXDIST;

  cost2=getDBT(EMPTY_TREE,reference_vertex);   
  for (i=1;i<=nj;i++)     
    {       
      int reference_child=T2->child(reference_vertex,i);
      dist_v_w_v_l = _d_v_w_v_l.getDBT(input_vertex,reference_child)-getDBT(EMPTY_TREE,reference_child);


      if (T2->childIsInComplex(reference_vertex,i))
	{
	  if (dist_v_w_v_l<min)
	    {
	      min = dist_v_w_v_l;
	      jm  = reference_child;
	      mtc = 2;
	    }
	}
    }               
  cost2=cost2+min; 

  if (cost2<MIN) 
    {
      MIN=cost2; 
      MTC=mtc; 
    }

  min = _d_v_w_v_l.getDBF(input_vertex,reference_vertex);
  cost3 = _d_l_l_v_l.getDBF(input_vertex,reference_vertex);
  if (min<=cost3)
    {
      mtc = 3;
      cost3 = min+_d_v_w_v_l.getCCost(input_vertex,reference_vertex);
    }
  else
    {
      mtc = 4;
      cost3 += _d_l_l_v_l.getCCost(input_vertex,reference_vertex);
    }

  if (cost3<MIN) 
    { 
      MIN=cost3;
      MTC=mtc; 
    }
  
  //-----------------------------------
  // We maintain the matching lists
  // mise a jour des listes d'alignement
  //-----------------------------------
  switch (MTC)
    {
    case 1 :{
      _choices_v_w_v_l.putFirst(input_vertex,reference_vertex,im); 
      _choices_v_w_v_l.putLast(input_vertex,reference_vertex,-1); 
    }
    break;
    case 2 :{
      _choices_v_w_v_l.putFirst(input_vertex,reference_vertex,jm);
      _choices_v_w_v_l.putLast(input_vertex,reference_vertex,M_v_w_v_l(input_vertex,jm));
    }
   
    break;

    case 3 : case 4 :{
      _choices_v_w_v_l.putFirst(input_vertex,reference_vertex,-1);
      _choices_v_w_v_l.putLast(input_vertex,reference_vertex,reference_vertex);
    }
    break;
    default: 	assert(0);break;
    }
  _choices_v_w_v_l.putFirst(input_vertex,reference_vertex,MTC);
  _d_v_w_v_l.putDBT(input_vertex,reference_vertex,MIN);
  return(MIN);
}

DistanceType MultiscaleMatching::distanceBetweenTree_v_w_l_w(int input_vertex,int reference_vertex)
{
  // On stocke dans ni et nj le nombre d'enfants de T1[i] et T2[j]
  int ni=T1->getNbChild(input_vertex);
  int nj=T2->getNbChild(reference_vertex);
  
  DistanceType cost1,cost2,cost3;
  DistanceType min,MIN=2*MAXDIST;
  int im=0,jm=0,MTC=0,mtc=0;
  int i;
  DistanceType dist,dist_v_w_l_w;

  min=MAXDIST;

  cost1=getDBT(input_vertex,EMPTY_TREE);   
  for (i=1;i<=ni;i++)
    {
      int input_child=T1->child(input_vertex,i);
      dist_v_w_l_w = _d_v_w_l_w.getDBT(input_child,reference_vertex)-getDBT(input_child,EMPTY_TREE);

      if (T1->childIsInComplex(input_vertex,i))
	{
	  if (dist_v_w_l_w<min)
	    {
	      min = dist_v_w_l_w;
	      im  = input_child;
	      mtc = 1;
	    }  
	}
    }
  cost1=cost1+min;
  // On conserve le cout minimum
  if (cost1<MIN)
    {
      MIN=cost1; 
      MTC=mtc; 
    }

  min=MAXDIST;

  cost2=getDBT(EMPTY_TREE,reference_vertex);   
  for (i=1;i<=nj;i++)     
    {       
      int reference_child=T2->child(reference_vertex,i);
      dist_v_w_l_w = _d_v_w_l_w.getDBT(input_vertex,reference_child)-getDBT(EMPTY_TREE,reference_child);


      if (T2->childIsInComplex(reference_vertex,i))
	{
	  if (dist_v_w_l_w<min)
	    {
	      min = dist_v_w_l_w;
	      jm  = reference_child;
	      mtc = 2;
	    }
	}
    }               
  cost2=cost2+min; 

  if (cost2<MIN) 
    {
      MIN=cost2; 
      MTC=mtc; 
    }

  min = _d_v_w_l_w.getDBF(input_vertex,reference_vertex);
  cost3 = _d_l_l_l_w.getDBF(input_vertex,reference_vertex);
  if (min<=cost3)
    {
      mtc = 3;
      cost3 = min+_d_v_w_l_w.getCCost(input_vertex,reference_vertex);
    }
  else
    {
      mtc = 4;
      cost3 += _d_l_l_l_w.getCCost(input_vertex,reference_vertex);
    }

  if (cost3<MIN) 
    { 
      MIN=cost3;
      MTC=mtc; 
    }
  
  //-----------------------------------
  // We maintain the matching lists
  // mise a jour des listes d'alignement
  //-----------------------------------
  switch (MTC)
    {
    case 1 :{
      _choices_v_w_l_w.putFirst(input_vertex,reference_vertex,im); 
      _choices_v_w_l_w.putLast(input_vertex,reference_vertex,-1); 
    }
    break;
    case 2 :{
      _choices_v_w_l_w.putFirst(input_vertex,reference_vertex,jm);
      _choices_v_w_l_w.putLast(input_vertex,reference_vertex,M_v_w_l_w(input_vertex,jm));
    }
   
    break;

    case 3 : case 4 :{
      _choices_v_w_l_w.putFirst(input_vertex,reference_vertex,-1);
      _choices_v_w_l_w.putLast(input_vertex,reference_vertex,reference_vertex);
    }
    break;
    default: 	assert(0);break;
    }
  _choices_v_w_l_w.putFirst(input_vertex,reference_vertex,MTC);
  _d_v_w_l_w.putDBT(input_vertex,reference_vertex,MIN);
  return(MIN);
}


DistanceType MultiscaleMatching::distanceBetweenTree_v_w_v_w(int input_vertex,int reference_vertex)
{
  // On stocke dans ni et nj le nombre d'enfants de T1[i] et T2[j]
  int ni=T1->getNbChild(input_vertex);
  int nj=T2->getNbChild(reference_vertex);
  
  DistanceType cost1,cost2,cost3;
  DistanceType min,MIN=2*MAXDIST;
  int im=0,jm=0,MTC=0,mtc=0;
  int i;
  DistanceType dist,dist_v_w_v_w;

  min=MAXDIST;

  cost1=getDBT(input_vertex,EMPTY_TREE);   
  for (i=1;i<=ni;i++)
    {
      int input_child=T1->child(input_vertex,i);
      dist_v_w_v_w = _d_v_w_v_w.getDBT(input_child,reference_vertex)-getDBT(input_child,EMPTY_TREE);

      if (T1->childIsInComplex(input_vertex,i))
	{
	  if (dist_v_w_v_w<min)
	    {
	      min = dist_v_w_v_w;
	      im  = input_child;
	      mtc = 1;
	    }  
	}
    }
  cost1=cost1+min;
  // On conserve le cout minimum
  if (cost1<MIN)
    {
      MIN=cost1; 
      MTC=mtc; 
    }

  min=MAXDIST;

  cost2=getDBT(EMPTY_TREE,reference_vertex);   
  for (i=1;i<=nj;i++)     
    {       
      int reference_child=T2->child(reference_vertex,i);
       dist_v_w_v_w = _d_v_w_v_w.getDBT(input_vertex,reference_child)-getDBT(EMPTY_TREE,reference_child);


      if (T2->childIsInComplex(reference_vertex,i))
	{
	  if (dist_v_w_v_w<min)
	    {
	      min = dist_v_w_v_w;
	      jm  = reference_child;
	      mtc = 2;
	    }
	}
    }               
  cost2=cost2+min; 

  if (cost2<MIN) 
    {
      MIN=cost2; 
      MTC=mtc; 
    }

  min = _d_v_w_v_w.getDBF(input_vertex,reference_vertex);
  cost3 = _d_l_l_v_w.getDBF(input_vertex,reference_vertex);
  if (min<=cost3)
    {
      mtc = 3;
      cost3 = min+_d_v_w_v_w.getCCost(input_vertex,reference_vertex);
    }
  else
    {
      mtc = 4;
      cost3 += _d_l_l_v_w.getCCost(input_vertex,reference_vertex);
    }

  if (cost3<MIN) 
    { 
      MIN=cost3;
      MTC=mtc; 
    }
  
  //-----------------------------------
  // We maintain the matching lists
  // mise a jour des listes d'alignement
  //-----------------------------------
  switch (MTC)
    {
    case 1 :{
      _choices_v_w_v_w.putFirst(input_vertex,reference_vertex,im); 
      _choices_v_w_v_w.putLast(input_vertex,reference_vertex,-1); 
    }
    break;
    case 2 :{
      _choices_v_w_v_w.putFirst(input_vertex,reference_vertex,jm);
      _choices_v_w_v_w.putLast(input_vertex,reference_vertex,M_v_w_v_w(input_vertex,jm));
    }
   
    break;

    case 3 : case 4 :{
      _choices_v_w_v_w.putFirst(input_vertex,reference_vertex,-1);
      _choices_v_w_v_w.putLast(input_vertex,reference_vertex,reference_vertex);
    }
    break;
    default: 	assert(0);break;
    }
  _choices_v_w_v_w.putFirst(input_vertex,reference_vertex,MTC);
  _d_v_w_v_w.putDBT(input_vertex,reference_vertex,MIN);
  return(MIN);
}


DistanceType MultiscaleMatching::distanceBetweenTree_v_w_(int input_vertex,int reference_vertex)
{
  // On stocke dans ni et nj le nombre d'enfants de T1[i] et T2[j]
  int ni=T1->getNbChild(input_vertex);
  int nj=T2->getNbChild(reference_vertex);
  
  DistanceType cost1,cost2,cost3;
  DistanceType min,MIN=2*MAXDIST;
  int im=0,jm=0,MTC=0,mtc=0;
  int i;
  DistanceType dist,dist_v_w_;

  min=MAXDIST;

  cost1=getDBT(input_vertex,EMPTY_TREE);   
  for (i=1;i<=ni;i++)
    {
      int input_child=T1->child(input_vertex,i);
      dist_v_w_ = _d_v_w_.getDBT(input_child,reference_vertex)-getDBT(input_child,EMPTY_TREE);

      if (T1->childIsInComplex(input_vertex,i))
	{
	  if (dist_v_w_<min)
	    {
	      min = dist_v_w_;
	      im  = input_child;
	      mtc = 1;
	    }  
	}
    }
  cost1=cost1+min;
  // On conserve le cout minimum
  if (cost1<MIN)
    {
      MIN=cost1; 
      MTC=mtc; 
    }

  min=MAXDIST;

  cost2=getDBT(EMPTY_TREE,reference_vertex);   
  for (i=1;i<=nj;i++)     
    {       
      int reference_child=T2->child(reference_vertex,i);
      dist_v_w_ = _d_v_w_.getDBT(input_vertex,reference_child)-getDBT(EMPTY_TREE,reference_child);


      if (T2->childIsInComplex(reference_vertex,i))
	{
	  if (dist_v_w_<min)
	    {
	      min = dist_v_w_;
	      jm  = reference_child;
	      mtc = 2;
	    }
	}
    }               
  cost2=cost2+min; 

  if (cost2<MIN) 
    {
      MIN=cost2; 
      MTC=mtc; 
    }

  cost3 = _d_v_w_.getDBF(input_vertex,reference_vertex);
  mtc = 3;
  cost3 += _d_v_w_.getCCost(input_vertex,reference_vertex);

  if (cost3<MIN) 
    { 
      MIN=cost3;
      MTC=mtc; 
    }
  
  //-----------------------------------
  // We maintain the matching lists
  // mise a jour des listes d'alignement
  //-----------------------------------
  switch (MTC)
    {
    case 1 :{
      _choices_v_w_.putFirst(input_vertex,reference_vertex,im); 
      _choices_v_w_.putLast(input_vertex,reference_vertex,-1); 
    }
    break;
    case 2 :{
      _choices_v_w_.putFirst(input_vertex,reference_vertex,jm);
      _choices_v_w_.putLast(input_vertex,reference_vertex,M_v_w_(input_vertex,jm));
    }
   
    break;

    case 3 : {
      _choices_v_w_.putFirst(input_vertex,reference_vertex,-1);
      _choices_v_w_.putLast(input_vertex,reference_vertex,reference_vertex);
    }
    break;
    default: 	assert(0);break;
    }
  _choices_v_w_.putFirst(input_vertex,reference_vertex,MTC);
  _d_v_w_.putDBT(input_vertex,reference_vertex,MIN);
  return(MIN);
}

DistanceType MultiscaleMatching::distanceBetweenTree_l_l_v_l(int input_vertex,int reference_vertex)
{
  
  DistanceType MIN;



  MIN = _d_l_l_v_l.getDBF(input_vertex,reference_vertex)+_d_l_l_v_l.getICost(reference_vertex)
    +_d_l_l_v_l.getDCost(input_vertex);
  
  _choices_l_l_v_l.putFirst(input_vertex,reference_vertex,-1);
  _choices_l_l_v_l.putLast(input_vertex,reference_vertex,-1);
  _choices_l_l_v_l.putFirst(input_vertex,reference_vertex,1);
  _d_l_l_v_l.putDBT(input_vertex,reference_vertex,MIN);
  return(MIN);
}
DistanceType MultiscaleMatching::distanceBetweenTree_l_l_l_w(int input_vertex,int reference_vertex)
{
   
  DistanceType MIN;



  MIN = _d_l_l_l_w.getDBF(input_vertex,reference_vertex)+_d_l_l_l_w.getICost(reference_vertex)
    +_d_l_l_l_w.getDCost(input_vertex);
  
  _choices_l_l_l_w.putFirst(input_vertex,reference_vertex,-1);
  _choices_l_l_l_w.putLast(input_vertex,reference_vertex,-1);
  _choices_l_l_l_w.putFirst(input_vertex,reference_vertex,1);
  _d_l_l_l_w.putDBT(input_vertex,reference_vertex,MIN);
  return(MIN);
}

DistanceType MultiscaleMatching::distanceBetweenTree_l_l_v_w(int input_vertex,int reference_vertex)
{
  
  DistanceType MIN;



  MIN = _d_l_l_v_w.getDBF(input_vertex,reference_vertex)+_d_l_l_v_w.getICost(reference_vertex)
    +_d_l_l_v_w.getDCost(input_vertex);
  
  _choices_l_l_v_w.putFirst(input_vertex,reference_vertex,-1);
  _choices_l_l_v_w.putLast(input_vertex,reference_vertex,-1);
  _choices_l_l_v_w.putFirst(input_vertex,reference_vertex,1);
  _d_l_l_v_w.putDBT(input_vertex,reference_vertex,MIN);
  return(MIN);
}

DistanceType MultiscaleMatching::distanceBetweenTree_l_l_l_l(int input_vertex,int reference_vertex)
{
  
  DistanceType MIN;

  MIN = _d_l_l_l_l.getDBF(input_vertex,reference_vertex)+_d_l_l_l_l.getICost(reference_vertex)
    +_d_l_l_l_l.getDCost(input_vertex);
  
  _choices_l_l_l_l.putFirst(input_vertex,reference_vertex,-1);
  _choices_l_l_l_l.putLast(input_vertex,reference_vertex,-1);
  _choices_l_l_l_l.putFirst(input_vertex,reference_vertex,1);
  _d_l_l_l_l.putDBT(input_vertex,reference_vertex,MIN);
  return(MIN);
}

DistanceType MultiscaleMatching::distanceBetweenTree_l_l(int input_vertex,int reference_vertex)
{
  
  DistanceType MIN=MAXDIST,dist_l_l_v_l,dist_l_l_l_w,dist_l_l_v_w,dist_l_l_l_l;
  int MTC =0;

  dist_l_l_v_l = distanceBetweenTree_l_l_v_l(input_vertex,reference_vertex);
  dist_l_l_l_w = distanceBetweenTree_l_l_l_w(input_vertex,reference_vertex);
  dist_l_l_v_w = distanceBetweenTree_l_l_v_w(input_vertex,reference_vertex);
  dist_l_l_l_l = distanceBetweenTree_l_l_l_l(input_vertex,reference_vertex);

  if (dist_l_l_v_l<MIN)
    {
      MIN = dist_l_l_v_l ;
      MTC = 1;
    }
  if (dist_l_l_l_w<MIN)
    {
      MIN = dist_l_l_l_w ;
      MTC = 2 ;
    }
  if (dist_l_l_v_w<MIN)
    {
      MIN = dist_l_l_v_w ;
      MTC =  3 ;
    }
  if (dist_l_l_l_l<MIN)
    {
      MIN = dist_l_l_l_l ;
      MTC =  4;
    }
   _choices_l_l.putFirst(input_vertex,reference_vertex,-1);
  _choices_l_l.putLast(input_vertex,reference_vertex,-1);
  _choices_l_l.putFirst(input_vertex,reference_vertex,MTC);
  _d_l_l.putDBT(input_vertex,reference_vertex,MIN);
  return(MIN);
}


DistanceType MultiscaleMatching::distanceBetweenTree_v_w(int input_vertex,int reference_vertex)
{
  
  DistanceType MIN=2*MAXDIST,dist_v_w_v_l,dist_v_w_l_w,dist_v_w_v_w,dist_v_w_;
  int MTC = 0, im=-1, jm= -1;

  dist_v_w_v_l = distanceBetweenTree_v_w_v_l(input_vertex,reference_vertex);
  dist_v_w_l_w = distanceBetweenTree_v_w_l_w(input_vertex,reference_vertex);
  dist_v_w_v_w = distanceBetweenTree_v_w_v_w(input_vertex,reference_vertex);
  dist_v_w_ = distanceBetweenTree_v_w_(input_vertex, reference_vertex);

  if (dist_v_w_v_l<MIN)
    {
      MIN = dist_v_w_v_l ;
      MTC = 1;
    }
  if (dist_v_w_l_w<MIN)
    {
      MIN = dist_v_w_l_w ;
      MTC = 2 ;
    }
  if (dist_v_w_v_w<MIN)
    {
      MIN = dist_v_w_v_w ;
      MTC =  3 ;
    }
  if (dist_v_w_<MIN)
    {
      MIN = dist_v_w_ ;
      MTC =  4;
    }
  _choices_v_w.putFirst(input_vertex,reference_vertex,-1); 
  _choices_v_w.putLast(input_vertex,reference_vertex,-1); 
  _choices_v_w.putFirst(input_vertex,reference_vertex,MTC);
  _d_v_w.putDBT(input_vertex,reference_vertex,MIN);
  return(MIN);

}

DistanceType MultiscaleMatching::distanceBetweenTree(int input_vertex,int reference_vertex)
{
  
  DistanceType MIN=2*MAXDIST,dist_v_l,dist_l_w,dist_v_w,dist_l_l;
  int MTC=0,im=-1,jm=-1;

  dist_v_l = distanceBetweenTree_v_l(input_vertex,reference_vertex);
  dist_l_w = distanceBetweenTree_l_w(input_vertex,reference_vertex);
  dist_v_w = distanceBetweenTree_v_w(input_vertex,reference_vertex);
  dist_l_l = distanceBetweenTree_l_l(input_vertex,reference_vertex);

  if (dist_v_l<MIN)
    {
      MIN = dist_v_l ;
      MTC = 1 ;
    }
  if (dist_l_w<MIN)
    {
      MIN = dist_l_w ;
      MTC =  2 ;
    }
  if (dist_v_w<MIN)
    {
      MIN = dist_v_w ;
      MTC =   3;
    }
  _choices.putFirst(input_vertex,reference_vertex,-1);
  _choices.putLast(input_vertex,reference_vertex,-1);

  _choices.putFirst(input_vertex,reference_vertex,MTC);
  _d.putDBT(input_vertex,reference_vertex,MIN);

  return(MIN);
}


