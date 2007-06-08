#include "multiscale_matching.h"
#include "dec_matrix.h"



DistanceType MultiscaleMatching::distanceBetweenForest_v_l(int input_vertex,int reference_vertex) 
{
  int ni=T1->getNbChild(input_vertex);
  int nj=T2->getNbChild(reference_vertex);
  DistanceType cost1,cost2,cost3;
  DistanceType min,DIST;
  int im=0,jm=0,im3=0,jm3=0,MFC=0,mfc=0;
  int i;
  DistanceType dist,dist_v_l,dist_v_w;
  
  DIST=MAXDIST;

  min=MAXDIST;
  cost1=getDBF(input_vertex,EMPTY_TREE); 
  for (i=1;i<=ni;i++)
    {
      int input_child=T1->child(input_vertex,i);
      dist_v_l = _d_v_l.getDBF(input_child,reference_vertex)-getDBF(input_child,EMPTY_TREE);

     if (!(T1->childIsInComplex(input_vertex,i)))
	{
	  if (dist_v_l<min)
	    {
	      min = dist_v_l;
	      im  = input_child;
	      mfc = 1;
	    }  
	}
    }
  cost1=cost1+min;
  // On conserve le cout minimum
  if (cost1<DIST)
    {
      DIST=cost1; 
      MFC=mfc; 
    }

  min=MAXDIST;
  cost2=getDBF(EMPTY_TREE,reference_vertex);
  for (i=1;i<=nj;i++)     
    {       
      int reference_child=T2->child(reference_vertex,i);
      dist_v_l=_d_v_l.getDBF(input_vertex,reference_child)-getDBF(EMPTY_TREE,reference_child);    
      dist_v_w=_d_v_w.getDBF(input_vertex,reference_child)-getDBF(EMPTY_TREE,reference_child);    
       if (dist_v_l<min) 
	{
	  min=dist_v_l;
	  jm=reference_child;
	  mfc = 2;
	}
       if (!(T2->childIsInComplex(reference_vertex,i)))
	{
	  if (dist_v_w<min)
	    {
	      min = dist_v_w;
	      jm  = reference_child;
	      mfc = 3;
	    }
	}
    }               
  cost2=cost2+min; 

  if (cost2<DIST) 
    {
      DIST=cost2; 
      MFC=mfc; 
    }
  
  cost3 = MAXDIST;
  if (ni==0)
    { 
      cost3=getDBF(EMPTY_TREE,reference_vertex);
      mfc = 4;
    }
  else
    {
      if (nj==0)
	{ 
	  cost3=getDBF(input_vertex,EMPTY_TREE); 
	  mfc = 4;
	}
      else
	{
	  DistanceType indelforest =getDBF(input_vertex,EMPTY_TREE)+getDBF(EMPTY_TREE,reference_vertex); 
	  cost3 = indelforest ;
	  for (int s1=1;s1<=ni;s1++) 
	    {
	      int input_child = T1->child(input_vertex,s1);
	      for (int s2=1;s2<=nj;s2++) 
		{
		  int reference_child = T2->child(reference_vertex,s2);
		  if (T1->childIsInComplex(input_vertex,s1))
		    {
		      dist_v_l =  _d_v_l.getDBT(input_child,reference_child)+indelforest
			-getDBT(input_child,EMPTY_TREE)-getDBT(EMPTY_TREE,reference_child);
		      dist_v_w =  _d_v_w.getDBT(input_child,reference_child)+indelforest
			-getDBT(input_child,EMPTY_TREE)-getDBT(EMPTY_TREE,reference_child);
		      if (dist_v_l<cost3)
			{
			  cost3 = dist_v_l ;
			  im3 = input_child;
			  jm3 = reference_child;
			  mfc = 4;
			}			      			  
		      if (!(T2->childIsInComplex(reference_vertex,s2)))
			{
			  if (dist_v_w<cost3)
			    {
			      cost3 = dist_v_w ;
			      im3 = input_child;
			      jm3 = reference_child;
			      mfc = 5;
			    }			      			  
			}
		    }
		}
	    }
	}
    }
  
  
  if (cost3<=DIST) 
    {
      DIST=cost3; 
      MFC=mfc;
    }
  
  _choices_v_l.createList(input_vertex,reference_vertex);
  _choices_v_l.putFirst(input_vertex,reference_vertex,MFC);
  
  
  switch(MFC)
    {
    case 1 :  
      {
	_choices_v_l.putLast(input_vertex,reference_vertex,im);
      }
      break;
    case 2 : case 3 :
      {
	_choices_v_l.putLast(input_vertex,reference_vertex,jm);
      }
      break;
    case 4 : case 5 : 

      {
	for (int i=1;i<=T1->getNbChild(input_vertex);i++)
	  {
	    int input_child = T1->child(input_vertex,i);
	    if (input_child==im3)
	      _choices_v_l.putLast(input_vertex,reference_vertex,jm3);
	    else
	      _choices_v_l.putLast(input_vertex,reference_vertex,EMPTY_NODE);
	  }
      }
      break;
    default : 	break;
    }
  _d_v_l.putDBF(input_vertex,reference_vertex,DIST);
  return(DIST);
}

DistanceType MultiscaleMatching::distanceBetweenForest_l_w(int input_vertex,int reference_vertex) 
{
  int ni=T1->getNbChild(input_vertex);
  int nj=T2->getNbChild(reference_vertex);
  DistanceType cost1,cost2,cost3;
  DistanceType min,DIST;
  int im=0,jm=0,im3=0,jm3=0,MFC=0,mfc=0;
  int i;
  DistanceType dist,dist_l_w,dist_v_w;
  
  DIST=MAXDIST;

  min=MAXDIST;
  cost1=getDBF(input_vertex,EMPTY_TREE); 
  for (i=1;i<=ni;i++)
    {
      int input_child=T1->child(input_vertex,i);
      dist_l_w = _d_l_w.getDBF(input_child,reference_vertex)-getDBF(input_child,EMPTY_TREE);
      dist_v_w = _d_v_w.getDBF(input_child,reference_vertex)-getDBF(input_child,EMPTY_TREE);

      if (dist_l_w<min)
	{
	  min = dist_l_w;
	  im  = input_child;
	  mfc = 1;
	}  

      if (!(T1->childIsInComplex(input_vertex,i)))
	{
	  if (dist_v_w<min)
	    {
	      min = dist_v_w;
	      im  = input_child;
	      mfc = 2;
	    }  
	}
    }
  cost1=cost1+min;
  // On conserve le cout minimum
  if (cost1<DIST)
    {
      DIST=cost1; 
      MFC=mfc; 
    }

  min=MAXDIST;
  cost2=getDBF(EMPTY_TREE,reference_vertex);
  for (i=1;i<=nj;i++)     
    {       
      int reference_child=T2->child(reference_vertex,i);
      dist_l_w=_d_l_w.getDBF(input_vertex,reference_child)-getDBF(EMPTY_TREE,reference_child);    
      if (!(T2->childIsInComplex(reference_vertex,i)))
	{
	  if (dist_l_w<min)
	    {
	      min = dist_l_w;
	      jm  = reference_child;
	      mfc = 3;
	    }
	}
    }               
  cost2=cost2+min; 

  if (cost2<DIST) 
    {
      DIST=cost2; 
      MFC=mfc; 
    }
  
  cost3 = MAXDIST;
  if (ni==0)
    { 
      cost3=getDBF(EMPTY_TREE,reference_vertex);
      mfc = 4;
    }
  else
    {
      if (nj==0)
	{ 
	  cost3=getDBF(input_vertex,EMPTY_TREE); 
	  mfc = 4;
	}
      else
	{
	  DistanceType indelforest =getDBF(input_vertex,EMPTY_TREE)+getDBF(EMPTY_TREE,reference_vertex); 
	  cost3 = indelforest ;
	  for (int s1=1;s1<=ni;s1++) 
	    {
	      int input_child = T1->child(input_vertex,s1);
	      for (int s2=1;s2<=nj;s2++) 
		{
		  int reference_child = T2->child(reference_vertex,s2);
		  if (T2->childIsInComplex(reference_vertex,s2))
		    {
		      dist_l_w =  _d_l_w.getDBT(input_child,reference_child)+indelforest
			-getDBT(input_child,EMPTY_TREE)-getDBT(EMPTY_TREE,reference_child);
		      dist_v_w =  _d_v_w.getDBT(input_child,reference_child)+indelforest
			-getDBT(input_child,EMPTY_TREE)-getDBT(EMPTY_TREE,reference_child);
		      if (dist_l_w<cost3)
			{
			  cost3 = dist_l_w ;
			  im3 = input_child;
			  jm3 = reference_child;
			  mfc = 4;
			}			      			  
		      if (!(T1->childIsInComplex(input_vertex,s1)))
			{
			  if (dist_v_w<cost3)
			    {
			      cost3 = dist_v_w ;
			      im3 = input_child;
			      jm3 = reference_child;
			      mfc = 5;
			    }			      			  
			}
		    }
		}
	    }
	}
    }
  
  
  if (cost3<=DIST) 
    {
      DIST=cost3; 
      MFC=mfc;
    }
  
  _choices_l_w.createList(input_vertex,reference_vertex);
  _choices_l_w.putFirst(input_vertex,reference_vertex,MFC);
  
  
  switch(MFC)
    {
    case 1 : case 2 :
      {
	_choices_l_w.putLast(input_vertex,reference_vertex,im);
      }
      break;
    case 3 :
      {
	_choices_l_w.putLast(input_vertex,reference_vertex,jm);
      }
      break;
    case 4 :     case 5 : 
      {
	for (int i=1;i<=T1->getNbChild(input_vertex);i++)
	  {
	    int input_child = T1->child(input_vertex,i);
	    if (input_child==im3)
	      _choices_l_w.putLast(input_vertex,reference_vertex,jm3);
	    else
	      _choices_l_w.putLast(input_vertex,reference_vertex,EMPTY_NODE);
	  }
      }
      break;
    default : 	break;
    }
  _d_l_w.putDBF(input_vertex,reference_vertex,DIST);
  return(DIST);
  
}


DistanceType MultiscaleMatching::distanceBetweenForest_v_w_v_l(int input_vertex,int reference_vertex) 
{
  int ni=T1->getNbChild(input_vertex);
  int nj=T2->getNbChild(reference_vertex);
  DistanceType cost1,cost2,cost3;
  DistanceType min,DIST;
  int im=0,jm=0,im3=0,jm3=0,MFC=0,mfc=0;
  int i;
  DistanceType dist,dist_v_w_v_l,dist_v_w_v_w;
  
  DIST=MAXDIST;

  min=MAXDIST;
  cost1=getDBF(input_vertex,EMPTY_TREE); 
  for (i=1;i<=ni;i++)
    {
      int input_child=T1->child(input_vertex,i);
      dist_v_w_v_l = _d_v_w_v_l.getDBF(input_child,reference_vertex)-getDBF(input_child,EMPTY_TREE);

     if (T1->childIsInComplex(input_vertex,i))
	{
	  if (dist_v_w_v_l<min)
	    {
	      min = dist_v_w_v_l;
	      im  = input_child;
	      mfc = 1;
	    }  
	}
    }
  cost1=cost1+min;
  // On conserve le cout minimum
  if (cost1<DIST)
    {
      DIST=cost1; 
      MFC=mfc; 
    }

  min=MAXDIST;
  cost2=getDBF(EMPTY_TREE,reference_vertex);
  for (i=1;i<=nj;i++)     
    {       
      int reference_child=T2->child(reference_vertex,i);
      dist_v_w_v_l=_d_v_w_v_l.getDBF(input_vertex,reference_child)-getDBF(EMPTY_TREE,reference_child);    

      if (T2->childIsInComplex(reference_vertex,i))
	{
	  if (dist_v_w_v_l<min) 
	    {
	      min=dist_v_w_v_l;
	      jm=reference_child;
	      mfc = 2;
	    }
	}
    }               
  cost2=cost2+min; 

  if (cost2<DIST) 
    {
      DIST=cost2; 
      MFC=mfc; 
    }
  
  cost3 = MAXDIST;
  if (ni==0)
    { 
      cost3=getDBF(EMPTY_TREE,reference_vertex);
      mfc = 3;
    }
  else
    {
      if (nj==0)
	{ 
	  cost3=getDBF(input_vertex,EMPTY_TREE); 
	  mfc = 3;
	}
      else
	{
	  DistanceType indelforest =getDBF(input_vertex,EMPTY_TREE)+getDBF(EMPTY_TREE,reference_vertex); 
	  cost3 = indelforest ;
	  for (int s1=1;s1<=ni;s1++) 
	    {
	      int input_child = T1->child(input_vertex,s1);
	      for (int s2=1;s2<=nj;s2++) 
		{
		  int reference_child = T2->child(reference_vertex,s2);
		  if ((T1->childIsInComplex(input_vertex,s1))&&(T2->childIsInComplex(reference_vertex,s2)))
		    {
		      dist_v_w_v_l =  _d_v_w_v_l.getDBT(input_child,reference_child)+indelforest
			-getDBT(input_child,EMPTY_TREE)-getDBT(EMPTY_TREE,reference_child);
		      if (dist_v_w_v_l<cost3)
			{
			  cost3 = dist_v_w_v_l ;
			  im3 = input_child;
			  jm3 = reference_child;
			  mfc = 3;
			}			      			  
		    }
		}
	    }
	}
    }
  
  
  if (cost3<=DIST) 
    {
      DIST=cost3; 
      MFC=mfc;
    }
  
  _choices_v_w_v_l.createList(input_vertex,reference_vertex);
  _choices_v_w_v_l.putFirst(input_vertex,reference_vertex,MFC);
  
  
  switch(MFC)
    {
    case 1 :  
      {
	_choices_v_w_v_l.putLast(input_vertex,reference_vertex,im);
      }
      break;
    case 2 : 
      {
	_choices_v_w_v_l.putLast(input_vertex,reference_vertex,jm);
      }
      break;
    case 3 : 
      {
	for (int i=1;i<=T1->getNbChild(input_vertex);i++)
	  {
	    int input_child = T1->child(input_vertex,i);
	    if (input_child==im3)
	      _choices_v_w_v_l.putLast(input_vertex,reference_vertex,jm3);
	    else
	      _choices_v_w_v_l.putLast(input_vertex,reference_vertex,EMPTY_NODE);
	  }
      }
      break;
    default : 	break;
    }
  _d_v_w_v_l.putDBF(input_vertex,reference_vertex,DIST);
  return(DIST);
}


DistanceType MultiscaleMatching::distanceBetweenForest_v_w_l_w(int input_vertex,int reference_vertex) 
{
  int ni=T1->getNbChild(input_vertex);
  int nj=T2->getNbChild(reference_vertex);
  DistanceType cost1,cost2,cost3;
  DistanceType min,DIST;
  int im=0,jm=0,im3=0,jm3=0,MFC=0,mfc=0;
  int i;
  DistanceType dist,dist_v_w_l_w,dist_v_w_v_w;
  
  DIST=MAXDIST;

  min=MAXDIST;
  cost1=getDBF(input_vertex,EMPTY_TREE); 
  for (i=1;i<=ni;i++)
    {
      int input_child=T1->child(input_vertex,i);
      dist_v_w_l_w = _d_v_w_l_w.getDBF(input_child,reference_vertex)-getDBF(input_child,EMPTY_TREE);

     if (T1->childIsInComplex(input_vertex,i))
	{
	  if (dist_v_w_l_w<min)
	    {
	      min = dist_v_w_l_w;
	      im  = input_child;
	      mfc = 1;
	    }  
	}
    }
  cost1=cost1+min;
  // On conserve le cout minimum
  if (cost1<DIST)
    {
      DIST=cost1; 
      MFC=mfc; 
    }

  min=MAXDIST;
  cost2=getDBF(EMPTY_TREE,reference_vertex);
  for (i=1;i<=nj;i++)     
    {       
      int reference_child=T2->child(reference_vertex,i);
      dist_v_w_l_w=_d_v_w_l_w.getDBF(input_vertex,reference_child)-getDBF(EMPTY_TREE,reference_child);    

      if (T2->childIsInComplex(reference_vertex,i))
	{
	  if (dist_v_w_l_w<min) 
	    {
	      min=dist_v_w_l_w;
	      jm=reference_child;
	      mfc = 2;
	    }
	}
    }               
  cost2=cost2+min; 

  if (cost2<DIST) 
    {
      DIST=cost2; 
      MFC=mfc; 
    }
  
  cost3 = MAXDIST;
  if (ni==0)
    { 
      cost3=getDBF(EMPTY_TREE,reference_vertex);
      mfc = 3;
    }
  else
    {
      if (nj==0)
	{ 
	  cost3=getDBF(input_vertex,EMPTY_TREE); 
	  mfc = 3;
	}
      else
	{
	  DistanceType indelforest =getDBF(input_vertex,EMPTY_TREE)+getDBF(EMPTY_TREE,reference_vertex); 
	  cost3 = indelforest ;
	  for (int s1=1;s1<=ni;s1++) 
	    {
	      int input_child = T1->child(input_vertex,s1);
	      for (int s2=1;s2<=nj;s2++) 
		{
		  int reference_child = T2->child(reference_vertex,s2);
		    if ((T1->childIsInComplex(input_vertex,s1))&&(T2->childIsInComplex(reference_vertex,s2)))
		    {
		      dist_v_w_l_w =  _d_v_w_l_w.getDBT(input_child,reference_child)+indelforest
			-getDBT(input_child,EMPTY_TREE)-getDBT(EMPTY_TREE,reference_child);
		      if (dist_v_w_l_w<cost3)
			{
			  cost3 = dist_v_w_l_w ;
			  im3 = input_child;
			  jm3 = reference_child;
			  mfc = 3;
			}			      			  
		    }
		}
	    }
	}
    }
  
  
  if (cost3<=DIST) 
    {
      DIST=cost3; 
      MFC=mfc;
    }
  
  _choices_v_w_l_w.createList(input_vertex,reference_vertex);
  _choices_v_w_l_w.putFirst(input_vertex,reference_vertex,MFC);
  
  
  switch(MFC)
    {
    case 1 :  
      {
	_choices_v_w_l_w.putLast(input_vertex,reference_vertex,im);
      }
      break;
    case 2 :
      {
	_choices_v_w_l_w.putLast(input_vertex,reference_vertex,jm);
      }
      break;
    case 3 : 
      {
	for (int i=1;i<=T1->getNbChild(input_vertex);i++)
	  {
	    int input_child = T1->child(input_vertex,i);
	    if (input_child==im3)
	      _choices_v_w_l_w.putLast(input_vertex,reference_vertex,jm3);
	    else
	      _choices_v_w_l_w.putLast(input_vertex,reference_vertex,EMPTY_NODE);
	  }
      }
      break;
    default : 	break;
    }
  _d_v_w_l_w.putDBF(input_vertex,reference_vertex,DIST);
  return(DIST);
}


DistanceType MultiscaleMatching::distanceBetweenForest_v_w_(int input_vertex,int reference_vertex) 
{
  int ni=T1->getNbChild(input_vertex);
  int nj=T2->getNbChild(reference_vertex);
  DistanceType cost1,cost2,cost3;
  DistanceType min,DIST;
  int im=0,jm=0,im3=0,jm3=0,MFC=0,mfc=0;
  int i;
  DistanceType dist,dist_v_w_,dist_v_w_v_w;
  
  DIST=MAXDIST;

  min=MAXDIST;
  cost1=getDBF(input_vertex,EMPTY_TREE); 
  for (i=1;i<=ni;i++)
    {
      int input_child=T1->child(input_vertex,i);
      dist_v_w_ = _d_v_w_.getDBF(input_child,reference_vertex)-getDBF(input_child,EMPTY_TREE);

     if (T1->childIsInComplex(input_vertex,i))
	{
	  if (dist_v_w_<min)
	    {
	      min = dist_v_w_;
	      im  = input_child;
	      mfc = 1;
	    }  
	}
    }
  cost1=cost1+min;
  // On conserve le cout minimum
  if (cost1<DIST)
    {
      DIST=cost1; 
      MFC=mfc; 
    }

  min=MAXDIST;
  cost2=getDBF(EMPTY_TREE,reference_vertex);
  for (i=1;i<=nj;i++)     
    {       
      int reference_child=T2->child(reference_vertex,i);
      dist_v_w_=_d_v_w_.getDBF(input_vertex,reference_child)-getDBF(EMPTY_TREE,reference_child);    

      if (T2->childIsInComplex(reference_vertex,i))
	{
	  if (dist_v_w_<min) 
	    {
	      min=dist_v_w_;
	      jm=reference_child;
	      mfc = 2;
	    }
	}
    }               
  cost2=cost2+min; 

  if (cost2<DIST) 
    {
      DIST=cost2; 
      MFC=mfc; 
    }
  
  cost3 = MAXDIST;
  NodeList* input_list=new NodeList();
  NodeList* reference_list=new NodeList();
  

  IntVector choice(ni+1);
 
  NEW_MAT(int,choice_v_w_,ni+1,nj+1)



  for (int s1=1;s1<=ni;s1++) 
    input_list->push_back(T1->child(input_vertex,s1));

  for (int s2=1;s2<=nj;s2++) 
    reference_list->push_back(T2->child(reference_vertex,s2)); 

 
  _restrMapp_v_w_.make(*input_list,*reference_list);
  _restrMappList_v_w_.resize(ni+nj+3,EMPTY_NODE);
   
  if (ni==0) 
    { 
      _restrMappList_v_w_[1]=2;
      for (i=1;i<=nj;i++)
	{
	  _restrMappList_v_w_[i+1]=1;
	}	  
      cost3=getDBF(EMPTY_TREE,reference_vertex);
   }

  else
    {
      if (nj==0) 
	{ 
	  _restrMappList_v_w_[2]=1;
	  for (i=1;i<=ni;i++)
	    _restrMappList_v_w_[i]=ni+1;	
	  cost3=getDBF(input_vertex,EMPTY_TREE); 
	}
      else
	{
	  DistanceType indelforest =getDBF(input_vertex,EMPTY_TREE)+getDBF(EMPTY_TREE,reference_vertex); 
	  cost3 = indelforest ;
	  for (int s1=1;s1<=ni;s1++) 
	    {
	      int input_child = T1->child(input_vertex,s1);
	      for (int s2=1;s2<=nj;s2++) 
		{
		  int reference_child = T2->child(reference_vertex,s2);
		  DistanceType indeltrees = getDBT(input_child,EMPTY_TREE)+getDBT(EMPTY_TREE,reference_child);
		  if ((T1->childIsInComplex(input_vertex,s1))&&(T2->childIsInComplex(reference_vertex,s2)))
		    {
		      dist_v_w_ =  _d_v_w_.getDBT(input_child,reference_child);
		      _restrDistances_v_w_.putDBT(input_child,reference_child,dist_v_w_);
		    }
		  else
		      _restrDistances_v_w_.putDBT(input_child,reference_child,indeltrees);
		}
	    }
	  cost3 = _restrMapp_v_w_.minCostFlow(_restrMappList_v_w_); 
	}
    }
  
  if (cost3<=DIST) 
    {
      DIST=cost3; 
      MFC=3  ;
    }
  
  _choices_v_w_.createList(input_vertex,reference_vertex);
  _choices_v_w_.putFirst(input_vertex,reference_vertex,MFC);
  
  
  switch(MFC)
    {
    case 1 :  
      {
	_choices_v_w_.putLast(input_vertex,reference_vertex,im);
      }
      break;
    case 2 :
      {
	_choices_v_w_.putLast(input_vertex,reference_vertex,jm);
      }
      break;
    case 3 : 
      {
	for (int i=1;i<=T1->getNbChild(input_vertex);i++)
	  {
	    int s = _restrMapp_v_w_.who(_restrMappList_v_w_[i]);
	    if (ni!=nj)
	      {
		if ((s != ni+1)&&(s != -1)) // alors s est un fils de reference_vertex
		  {
		    int input_child = T1->child(input_vertex,i);
		    int reference_child = s;
		    DistanceType indeltrees = getDBT(input_child,EMPTY_TREE)+getDBT(EMPTY_TREE,reference_child);	    
		    if (_restrDistances_v_w_.getDBT(input_child,reference_child)==indeltrees)
		      _choices_v_w_.putLast(input_vertex,reference_vertex,ni+1);
		    else
		      _choices_v_w_.putLast(input_vertex,reference_vertex,s);
		  } 
		else
		  _choices_v_w_.putLast(input_vertex,reference_vertex,s);
	      }
	    else 
	      {
		if (s != -1)
		  {
		    int input_child = T1->child(input_vertex,i);
		    int reference_child = s;
		    DistanceType indeltrees = getDBT(input_child,EMPTY_TREE)+getDBT(EMPTY_TREE,reference_child);	    
		    if (_restrDistances_v_w_.getDBT(input_child,reference_child)==indeltrees)
		      _choices_v_w_.putLast(input_vertex,reference_vertex,-1);
		    else
		      _choices_v_w_.putLast(input_vertex,reference_vertex,s);
		  }
		else
		  _choices_v_w_.putLast(input_vertex,reference_vertex,s);
	      } 
	  }
      }
      break;
    default : 	break;
    }
  _d_v_w_.putDBF(input_vertex,reference_vertex,DIST);
  DEL_MAT(choice_v_w_,ni+1)
  return(DIST);
}


DistanceType MultiscaleMatching::distanceBetweenForest_v_w_v_w(int input_vertex,int reference_vertex) 
{
  int ni=T1->getNbChild(input_vertex);
  int nj=T2->getNbChild(reference_vertex);
  DistanceType cost1,cost2,cost3;
  DistanceType min,DIST;
  int im=0,jm=0,im3=0,jm3=0,MFC=0,mfc=0;
  int i,j;
  DistanceType dist,dist_v_w_v_w,dist_v_w_l_w,dist_v_w_v_l,dist_v_l,dist_l_w,dist_v_w_;
  DistanceType dist_l_l_v_w,dist_l_l_l_l,dist_l_l;
  DIST=MAXDIST;

  min=MAXDIST;
  cost1=getDBF(input_vertex,EMPTY_TREE); 
  for (i=1;i<=ni;i++)
    {
      int input_child=T1->child(input_vertex,i);
      dist_v_w_v_w = _d_v_w_v_w.getDBF(input_child,reference_vertex)-getDBF(input_child,EMPTY_TREE);

     if (T1->childIsInComplex(input_vertex,i))
	{
	  if (dist_v_w_v_w<min)
	    {
	      min = dist_v_w_v_w;
	      im  = input_child;
	      mfc = 1;
	    }  
	}
    }
  cost1=cost1+min;
  // On conserve le cout minimum
  if (cost1<DIST)
    {
      DIST=cost1; 
      MFC=mfc; 
    }

  min=MAXDIST;
  cost2=getDBF(EMPTY_TREE,reference_vertex);
  for (i=1;i<=nj;i++)     
    {       
      int reference_child=T2->child(reference_vertex,i);
      dist_v_w_v_w=_d_v_w_v_w.getDBF(input_vertex,reference_child)-getDBF(EMPTY_TREE,reference_child);    

      if (T2->childIsInComplex(reference_vertex,i))
	{
	  if (dist_v_w_v_w<min) 
	    {
	      min=dist_v_w_v_w;
	      jm=reference_child;
	      mfc = 2;
	    }
	}
    }               
  cost2=cost2+min; 

  if (cost2<DIST) 
    {
      DIST=cost2; 
      MFC=mfc; 
    }
  
  cost3 = MAXDIST;
  NodeList* input_list=new NodeList();
  NodeList* reference_list=new NodeList();
  NEW_MAT(int,choice,ni+1,nj+1)
  IntVector choice_v_w_v_w;
  choice_v_w_v_w.resize(ni+1);
  
  for (int s1=1;s1<=ni;s1++) 
    input_list->push_back(T1->child(input_vertex,s1));

  for (int s2=1;s2<=nj;s2++) 
    reference_list->push_back(T2->child(reference_vertex,s2)); 

 
  _restrMapp_v_w_v_w.make(*input_list,*reference_list);
  _restrMappList_v_w_v_w.resize(ni+nj+3,EMPTY_NODE);
   
  if (ni==0) 
    { 
      _restrMappList_v_w_v_w[1]=2;
      for (i=1;i<=nj;i++)
	{
	  _restrMappList_v_w_v_w[i+1]=1;
	}	  
      cost3=getDBF(EMPTY_TREE,reference_vertex);
   }

  else
    {
      if (nj==0) 
	{ 
	  _restrMappList_v_w_v_w[2]=1;
	  for (i=1;i<=ni;i++)
	    _restrMappList_v_w_v_w[i]=ni+1;	
	  cost3=getDBF(input_vertex,EMPTY_TREE); 
	}
      else
	{
	  DistanceType indelforest =getDBF(input_vertex,EMPTY_TREE)+getDBF(EMPTY_TREE,reference_vertex); 
	  cost3 = indelforest ;
	  int s1;
	  for (s1=1;s1<=ni;s1++) 
	    {
	      int input_child = T1->child(input_vertex,s1);
	      for (int s2=1;s2<=nj;s2++) 
		{
		  int reference_child = T2->child(reference_vertex,s2);
		  dist = MAXDIST;
		  dist_v_l =  _d_v_l.getDBT(input_child,reference_child);
		  dist_l_w =  _d_l_w.getDBT(input_child,reference_child);
		  dist_l_l =  _d_l_l.getDBT(input_child,reference_child);
		  dist_v_w_ =  _d_v_w_.getDBT(input_child,reference_child);
		  dist_v_w_v_w =  _d_v_w_v_w.getDBT(input_child,reference_child);
		  dist_l_l_v_w =  _d_l_l_v_w.getDBT(input_child,reference_child);
		  dist_l_l_l_l =  _d_l_l_l_l.getDBT(input_child,reference_child);
		  if  (T1->childIsInComplex(input_vertex,s1))
		    {
		      if (T2->childIsInComplex(reference_vertex,s2))
			{
			  if (dist_v_w_v_w<dist)
			    {
			      choice[s1][s2]=1;
			      dist = dist_v_w_v_w;
			    }
			  if (dist_v_w_<dist)
			    {
			      choice[s1][s2]=2;
			      dist = dist_v_w_;
			    }
			  if (dist_l_l_v_w<dist)
			    {
			      choice[s1][s2]=3;
			      dist = dist_l_l_v_w;
			    }
			  if (dist_l_l_l_l<dist)
			    {
			      choice[s1][s2]=4;
			      dist = dist_l_l_l_l;
			    }
			  _restrDistances_v_w_v_w.putDBT(input_child,reference_child,dist);
			}
		      else
			{
			  if (dist_l_w<dist)
			    {
			      choice[s1][s2]=5;
			      dist = dist_l_w;
			    }
			  if (dist_l_l<dist)
			    {
			      choice[s1][s2]=7;
			      dist = dist_l_l;
			    }
			  _restrDistances_v_w_v_w.putDBT(input_child,reference_child,dist);
			}
		    }
		  else
		    {
		      if (T2->childIsInComplex(reference_vertex,s2))
			{
			  if (dist_v_l<dist)
			    {
			      choice[s1][s2]=6;
			      dist = dist_v_l;
			    }
			  if (dist_l_l<dist)
			    {
			      choice[s1][s2]=7;
			      dist = dist_l_l;
			    }
			  _restrDistances_v_w_v_w.putDBT(input_child,reference_child,dist);
			}
		      else
			{
			  dist =  _d.getDBT(input_child,reference_child);
			  choice[s1][s2]=8;
			  _restrDistances_v_w_v_w.putDBT(input_child,reference_child,dist);
			}
		    }
		}
	    }
	  DistanceVector temp;
	  temp.resize(nj+1);
	  for ( s1=1;s1<=ni;s1++) 
	    {
	      // _restrDistances = _restrDistances_v_w_v_w;
	      int input_child = T1->child(input_vertex,s1);
	      // Modification du graphe de flot
	      if (T1->childIsInComplex(input_vertex,s1))
		{
		  for (int s2=1;s2<=nj;s2++) 
		    {
		      int reference_child = T2->child(reference_vertex,s2);
		      if (!T2->childIsInComplex(reference_vertex,s2))
			{
			  temp[s2]=_restrDistances_v_w_v_w.getDBT(input_child,reference_child);
			  DistanceType indeltrees = getDBT(input_child,EMPTY_TREE)+getDBT(EMPTY_TREE,reference_child);
			  _restrDistances_v_w_v_w.putDBT(input_child,reference_child,indeltrees);
			  
			}
		    }
		  // calcul du flot optimal		  
		  min = _restrMapp_v_w_v_w.minCostFlow(_restrMappList_v_w_v_w);
		  if (min<cost3)
		    {
		      cost3 = min;
		      for (i=1;i<=ni;i++)
			{
			  choice_v_w_v_w[i] = _restrMapp_v_w_v_w.who(_restrMappList_v_w_v_w[i]);
			  int input_child = T1->child(input_vertex,i);
			  if (ni!=nj)
			    {
			      if ((choice_v_w_v_w[i]!=ni+1)&&(choice_v_w_v_w[i] != -1))
				{
				  int reference_child = choice_v_w_v_w[i];
				  DistanceType indeltrees = getDBT(input_child,EMPTY_TREE)+getDBT(EMPTY_TREE,reference_child);
				  if (_restrDistances_v_w_.getDBT(input_child,reference_child)==indeltrees)
				    choice_v_w_v_w[i] = ni+1 ;
				}
			    }
			  else
			    {
			      if (choice_v_w_v_w[i] != -1)
				{
				  int reference_child = choice_v_w_v_w[i];
				  DistanceType indeltrees = getDBT(input_child,EMPTY_TREE)+getDBT(EMPTY_TREE,reference_child);
				  if (_restrDistances_v_w_.getDBT(input_child,reference_child)==indeltrees)
				    choice_v_w_v_w[i] = -1 ;
				}
			    }
			}
		    }
		  // remise à jour du graphe de flot 
		  for (int s2b=1;s2b<=nj;s2b++) 
		    {
		      int reference_child = T2->child(reference_vertex,s2b);
		      if (!T2->childIsInComplex(reference_vertex,s2b))
			{
			  _restrDistances_v_w_v_w.putDBT(input_child,reference_child,temp[s2b]);
			  
			}
		    }
		}
	    }
	}
    }

  
  if (cost3<=DIST) 
    {
      DIST=cost3; 
      MFC=3;
    }
  
  _choices_v_w_v_w.createList(input_vertex,reference_vertex);
  _choices_v_w_v_w.putFirst(input_vertex,reference_vertex,MFC);
  
  
  switch(MFC)
    {
    case 1 :  
      {
	_choices_v_w_v_w.putLast(input_vertex,reference_vertex,im);
      }
      break;
    case 2 : 
      {
	_choices_v_w_v_w.putLast(input_vertex,reference_vertex,jm);
      }
      break;
    case 3: 
      {
	for (int i=1;i<=T1->getNbChild(input_vertex);i++)
	  {
	    if (choice_v_w_v_w[i]==-1)
	      {
		_choices_v_w_v_w.putLast(input_vertex,reference_vertex,-1);
	      }
	    else
	      {
		for (int j=1;j<=T2->getNbChild(reference_vertex);j++)
		  {
		    if (T2->child(reference_vertex,j)==choice_v_w_v_w[i])
		      {
			_choices_v_w_v_w.putLast(input_vertex,reference_vertex,choice[i][j]);
		      }
		  }
	      }
	    _choices_v_w_v_w.putLast(input_vertex,reference_vertex,choice_v_w_v_w[i]);
	  }
      }
      break;
    default : 	break;
    }
  _d_v_w_v_w.putDBF(input_vertex,reference_vertex,DIST);
  DEL_MAT(choice,ni+1)
  return(DIST);
}


DistanceType MultiscaleMatching::distanceBetweenForest_l_l_v_l(int input_vertex,int reference_vertex) 
{
  int ni=T1->getNbChild(input_vertex);
  int nj=T2->getNbChild(reference_vertex);
  DistanceType cost1,cost2,cost3;
  DistanceType min,DIST;
  int im=0,jm=0,im3=0,jm3=0,MFC=0,mfc=0;
  int i,j;
  DistanceType dist,dist_l_l_v_l,dist_l_l_v_w;
  
  DIST=MAXDIST;

  min=MAXDIST;
  cost1=getDBF(input_vertex,EMPTY_TREE); 
  for (i=1;i<=ni;i++)
    {
      int input_child=T1->child(input_vertex,i);
      dist_l_l_v_l = _d_l_l_v_l.getDBF(input_child,reference_vertex)-getDBF(input_child,EMPTY_TREE);

      if (!(T1->childIsInComplex(input_vertex,i)))
	{
	  if (dist_l_l_v_l<min)
	    {
	      min = dist_l_l_v_l;
	      im  = input_child;
	      mfc = 1;
	    }  
	}
    }
  cost1=cost1+min;
  // On conserve le cout minimum
  if (cost1<DIST)
    {
      DIST=cost1; 
      MFC=mfc; 
    }

  min=MAXDIST;
  cost2=getDBF(EMPTY_TREE,reference_vertex);
  for (i=1;i<=nj;i++)     
    {       
      int reference_child=T2->child(reference_vertex,i);
      dist_l_l_v_l=_d_l_l_v_l.getDBF(input_vertex,reference_child)-getDBF(EMPTY_TREE,reference_child);    
      dist_l_l_v_w=_d_l_l_v_w.getDBF(input_vertex,reference_child)-getDBF(EMPTY_TREE,reference_child);    

      if (dist_l_l_v_l<min) 
	{
	  min=dist_l_l_v_l;
	  jm=reference_child;
	  mfc = 2;
	}
      if (!(T2->childIsInComplex(reference_vertex,i)))
	{
	  if (dist_l_l_v_w<min) 
	    {
	      min=dist_l_l_v_w;
	      jm=reference_child;
	      mfc = 3;
	    }
	}
    }               
  cost2=cost2+min; 

  if (cost2<DIST) 
    {
      DIST=cost2; 
      MFC=mfc; 
    }
  
  cost3 = MAXDIST;
  if (ni==0)
    { 
      cost3=getDBF(EMPTY_TREE,reference_vertex);
      mfc = 4;
    }
  else
    {
      if (nj==0)
	{ 
	  cost3=getDBF(input_vertex,EMPTY_TREE); 
	  mfc = 4;
	}
      else
	{
	  DistanceType indelforest =getDBF(input_vertex,EMPTY_TREE)+getDBF(EMPTY_TREE,reference_vertex); 
	  cost3 = indelforest ;
	  for (int s1=1;s1<=ni;s1++) 
	    {
	      int input_child = T1->child(input_vertex,s1);
	      for (int s2=1;s2<=nj;s2++) 
		{
		  int reference_child = T2->child(reference_vertex,s2);
		  if (T1->childIsInComplex(input_vertex,s1))
		    {
		      dist_l_l_v_l =  _d_l_l_v_l.getDBT(input_child,reference_child)+indelforest
			-getDBT(input_child,EMPTY_TREE)-getDBT(EMPTY_TREE,reference_child);
		      if (dist_l_l_v_l<cost3)
			{
			  cost3 = dist_l_l_v_l ;
			  im3 = input_child;
			  jm3 = reference_child;
			  mfc = 4;
			}
		      if (!(T2->childIsInComplex(reference_vertex,s2)))
			{
			  dist_l_l_v_w =  _d_l_l_v_w.getDBT(input_child,reference_child)+indelforest
			    -getDBT(input_child,EMPTY_TREE)-getDBT(EMPTY_TREE,reference_child);
			  if (dist_l_l_v_w<cost3)
			    {
			      cost3 = dist_l_l_v_w ;
			      im3 = input_child;
			      jm3 = reference_child;
			      mfc = 5;
			    }
			}
		    }
		}
	    }
	}
    }
  
  
  if (cost3<=DIST) 
    {
      DIST=cost3; 
      MFC=mfc;
    }
  
  _choices_l_l_v_l.createList(input_vertex,reference_vertex);
  _choices_l_l_v_l.putFirst(input_vertex,reference_vertex,MFC);
  
  
  switch(MFC)
    {
    case 1 :  
      {
	_choices_l_l_v_l.putLast(input_vertex,reference_vertex,im);
      }
      break;
    case 2 : case 3 :
      {
	_choices_l_l_v_l.putLast(input_vertex,reference_vertex,jm);
      }
      break;
    case 4 : case 5 :
      {
	for (int i=1;i<=T1->getNbChild(input_vertex);i++)
	  {
	    int input_child = T1->child(input_vertex,i);
	    if (input_child==im3)
	      _choices_l_l_v_l.putLast(input_vertex,reference_vertex,jm3);
	    else
	      _choices_l_l_v_l.putLast(input_vertex,reference_vertex,EMPTY_NODE);
	  }
      }
      break;
    default : 	break;
    }
  _d_l_l_v_l.putDBF(input_vertex,reference_vertex,DIST);
  return(DIST);
}


DistanceType MultiscaleMatching::distanceBetweenForest_l_l_l_w(int input_vertex,int reference_vertex) 
{
  int ni=T1->getNbChild(input_vertex);
  int nj=T2->getNbChild(reference_vertex);
  DistanceType cost1,cost2,cost3;
  DistanceType min,DIST;
  int im=0,jm=0,im3=0,jm3=0,MFC=0,mfc=0;
  int i;
  DistanceType dist,dist_l_l_l_w,dist_l_l_l_l,dist_l_l_v_w;
  
  DIST=MAXDIST;

  min=MAXDIST;
  cost1=getDBF(input_vertex,EMPTY_TREE); 
  for (i=1;i<=ni;i++)
    {
      int input_child=T1->child(input_vertex,i);
      dist_l_l_l_w = _d_l_l_l_w.getDBF(input_child,reference_vertex)-getDBF(input_child,EMPTY_TREE);
      dist_l_l_v_w = _d_l_l_v_w.getDBF(input_child,reference_vertex)-getDBF(input_child,EMPTY_TREE);

      if (dist_l_l_l_w<min)
	{
	  min = dist_l_l_l_w;
	  im  = input_child;
	  mfc = 1;
	}  

      if (!(T1->childIsInComplex(input_vertex,i)))
	{
	  if (dist_l_l_v_w<min)
	    {
	      min = dist_l_l_v_w;
	      im  = input_child;
	      mfc = 2;
	    }  
	}
    }
  cost1=cost1+min;
  // On conserve le cout minimum
  if (cost1<DIST)
    {
      DIST=cost1; 
      MFC=mfc; 
    }

  min=MAXDIST;
  cost2=getDBF(EMPTY_TREE,reference_vertex);
  for (i=1;i<=nj;i++)     
    {       
      int reference_child=T2->child(reference_vertex,i);
      dist_l_l_l_w=_d_l_l_l_w.getDBF(input_vertex,reference_child)-getDBF(EMPTY_TREE,reference_child);    

      if (!(T2->childIsInComplex(reference_vertex,i)))
	{
	  if (dist_l_l_l_w<min) 
	    {
	      min=dist_l_l_l_w;
	      jm=reference_child;
	      mfc = 3;
	    }
	}
    }               
  cost2=cost2+min; 

  if (cost2<DIST) 
    {
      DIST=cost2; 
      MFC=mfc; 
    }
  
  cost3 = MAXDIST;
  if (ni==0)
    { 
      cost3=getDBF(EMPTY_TREE,reference_vertex);
      mfc = 4;
    }
  else
    {
      if (nj==0)
	{ 
	  cost3=getDBF(input_vertex,EMPTY_TREE); 
	  mfc = 4;
	}
      else
	{
	  DistanceType indelforest =getDBF(input_vertex,EMPTY_TREE)+getDBF(EMPTY_TREE,reference_vertex); 
	  cost3 = indelforest ;
	  for (int s1=1;s1<=ni;s1++) 
	    {
	      int input_child = T1->child(input_vertex,s1);
	      for (int s2=1;s2<=nj;s2++) 
		{
		  int reference_child = T2->child(reference_vertex,s2);
		    if (T2->childIsInComplex(reference_vertex,s2))
		    {
		      dist_l_l_l_w =  _d_l_l_l_w.getDBT(input_child,reference_child)+indelforest
			-getDBT(input_child,EMPTY_TREE)-getDBT(EMPTY_TREE,reference_child);
		      if (dist_l_l_l_w<cost3)
			{
			  cost3 = dist_l_l_l_w ;
			  im3 = input_child;
			  jm3 = reference_child;
			  mfc = 4;
			}		
		      if (!(T1->childIsInComplex(input_vertex,s1)))
			{
			  dist_l_l_v_w =  _d_l_l_v_w.getDBT(input_child,reference_child)+indelforest
			    -getDBT(input_child,EMPTY_TREE)-getDBT(EMPTY_TREE,reference_child);
			  if (dist_l_l_v_w<cost3)
			    {
			      cost3 = dist_l_l_v_w ;
			      im3 = input_child;
			      jm3 = reference_child;
			      mfc = 5;
			    }		
			}	      			  
		    }
		}
	    }
	}
    }
  
  
  if (cost3<=DIST) 
    {
      DIST=cost3; 
      MFC=mfc;
    }
  
  _choices_l_l_l_w.createList(input_vertex,reference_vertex);
  _choices_l_l_l_w.putFirst(input_vertex,reference_vertex,MFC);
  
  
  switch(MFC)
    {
    case 1 : case 2 : 
      {
	_choices_l_l_l_w.putLast(input_vertex,reference_vertex,im);
      }
      break;
   case 3 :
      {
	_choices_l_l_l_w.putLast(input_vertex,reference_vertex,jm);
      }
      break;
    case 4 :  case 5 : 
      {
	for (int i=1;i<=T1->getNbChild(input_vertex);i++)
	  {
	    int input_child = T1->child(input_vertex,i);
	    if (input_child==im3)
	      _choices_l_l_l_w.putLast(input_vertex,reference_vertex,jm3);
	    else
	      _choices_l_l_l_w.putLast(input_vertex,reference_vertex,EMPTY_NODE);
	  }
      }
      break;
    default : 	break;
    }
  _d_l_l_l_w.putDBF(input_vertex,reference_vertex,DIST);
  return(DIST);
}


DistanceType MultiscaleMatching::distanceBetweenForest_l_l_l_l(int input_vertex,int reference_vertex) 
{
  int ni=T1->getNbChild(input_vertex);
  int nj=T2->getNbChild(reference_vertex);
  DistanceType cost1,cost2,cost3;
  DistanceType min,DIST;
  int im=0,jm=0,im3=0,jm3=0,MFC=0,mfc=0;
  int i;
  DistanceType dist,dist_l_l_l_l,dist_l_l_v_l,dist_l_l_l_w,dist_l_l,dist_l_w,dist_v_l;
  
  DIST=MAXDIST;

  min=MAXDIST;
  cost1=getDBF(input_vertex,EMPTY_TREE); 
  for (i=1;i<=ni;i++)
    {
      int input_child=T1->child(input_vertex,i);
      dist_l_l_l_l = _d_l_l_l_l.getDBF(input_child,reference_vertex)-getDBF(input_child,EMPTY_TREE);
      dist_l_l_v_l = _d_l_l_v_l.getDBF(input_child,reference_vertex)-getDBF(input_child,EMPTY_TREE);

      if (dist_l_l_l_l<min)
	{
	  min = dist_l_l_l_l;
	  im  = input_child;
	  mfc = 1;
	}  
      if (!(T1->childIsInComplex(input_vertex,i)))
	{
	  if (dist_l_l_v_l<min)
	    {
	      min = dist_l_l_v_l;
	      im  = input_child;
	      mfc = 2;
	    }  
	}
    }
  cost1=cost1+min;
  // On conserve le cout minimum
  if (cost1<DIST)
    {
      DIST=cost1; 
      MFC=mfc; 
    }

  min=MAXDIST;
  cost2=getDBF(EMPTY_TREE,reference_vertex);
  for (i=1;i<=nj;i++)     
    {       
      int reference_child=T2->child(reference_vertex,i);
      dist_l_l_l_l=_d_l_l_l_l.getDBF(input_vertex,reference_child)-getDBF(EMPTY_TREE,reference_child);    
      dist_l_l_l_w=_d_l_l_l_w.getDBF(input_vertex,reference_child)-getDBF(EMPTY_TREE,reference_child);    

      if (dist_l_l_l_l<min) 
	{
	  min=dist_l_l_l_l;
	  jm=reference_child;
	  mfc = 3;
	}
      if (T2->childIsInComplex(reference_vertex,i))
	{
	  if (dist_l_l_l_w<min) 
	    {
	      min=dist_l_l_l_w;
	      jm=reference_child;
	      mfc = 4;
	    }
	}
    }               
  cost2=cost2+min; 

  if (cost2<DIST) 
    {
      DIST=cost2; 
      MFC=mfc; 
    }
  
  cost3 = MAXDIST;
  if (ni==0)
    { 
      cost3=getDBF(EMPTY_TREE,reference_vertex);
      mfc = 5;
    }
  else
    {
      if (nj==0)
	{ 
	  cost3=getDBF(input_vertex,EMPTY_TREE); 
	  mfc = 5;
	}
      else
	{
	  DistanceType indelforest =getDBF(input_vertex,EMPTY_TREE)+getDBF(EMPTY_TREE,reference_vertex); 
	  cost3 = indelforest ;
	  for (int s1=1;s1<=ni;s1++) 
	    {
	      int input_child = T1->child(input_vertex,s1);
	      for (int s2=1;s2<=nj;s2++) 
		{
		  int reference_child = T2->child(reference_vertex,s2);
		  dist_l_l_l_l =  _d_l_l_l_l.getDBT(input_child,reference_child)+indelforest
		    -getDBT(input_child,EMPTY_TREE)-getDBT(EMPTY_TREE,reference_child);
		  if (T1->childIsInComplex(input_vertex,s1))
		    {
		      if (T2->childIsInComplex(reference_vertex,s2))
			{
			  if (dist_l_l_l_l<cost3)
			    {
			      cost3 = dist_l_l_l_l ;
			      im3 = input_child;
			      jm3 = reference_child;
			      mfc = 5;
			    }	
			}
		      else // reference child is not in Complex
			{
			  dist_l_l_v_l =  _d_l_l_l_w.getDBT(input_child,reference_child)+indelforest
			    -getDBT(input_child,EMPTY_TREE)-getDBF(EMPTY_TREE,reference_child);
			  dist_l_w =  _d_l_w.getDBT(input_child,reference_child)+indelforest
			    -getDBT(input_child,EMPTY_TREE)-getDBT(EMPTY_TREE,reference_child);
			  if (dist_l_w<cost3)
			    {
			      cost3 = dist_l_w ;
			      im3 = input_child;
			      jm3 = reference_child;
			      mfc = 6;
			    }		
			  if (dist_l_l_l_l<cost3)
			    {
			      cost3 = dist_l_l_l_l ;
			      im3 = input_child;
			      jm3 = reference_child;
			      mfc = 5;
			    }
			  if (dist_l_l_v_l<cost3)
			    {
			      cost3 = dist_l_l_v_l ;
			      im3 = input_child;
			      jm3 = reference_child;
			      mfc = 8;
			    }	
			}
		    }
		  else // input  child is not in Complex
		    {
		      if (T2->childIsInComplex(reference_vertex,s2))
			{
			  dist_v_l =  _d_v_l.getDBT(input_child,reference_child)+indelforest
			    -getDBT(input_child,EMPTY_TREE)-getDBF(EMPTY_TREE,reference_child);
			  dist_l_l_l_w =  _d_l_l_l_w.getDBT(input_child,reference_child)+indelforest
			    -getDBT(input_child,EMPTY_TREE)-getDBT(EMPTY_TREE,reference_child);
			  if (dist_v_l<cost3)
			    {
			      cost3 = dist_v_l ;
			      im3 = input_child;
			      jm3 = reference_child;
			      mfc = 7;
			    }
			  if (dist_l_l_l_l<cost3)
			    {
			      cost3 = dist_l_l_l_l ;
			      im3 = input_child;
			      jm3 = reference_child;
			      mfc = 5;
			    }
			  if (dist_l_l_l_w<cost3)
			    {
			      cost3 = dist_l_l_l_w ;
			      im3 = input_child;
			      jm3 = reference_child;
			      mfc = 9;
			    }	
			}
		      else
			{
			  dist =  _d.getDBT(input_child,reference_child)+indelforest
			    -getDBT(input_child,EMPTY_TREE)-getDBT(EMPTY_TREE,reference_child);
			  if (dist<cost3)
			    {
			      cost3 = dist ;
			      im3 = input_child;
			      jm3 = reference_child;
			      mfc = 10;
			    }	
			}		        
			
		    }
		}
	    }
	}
    }
    
  if (cost3<=DIST) 
    {
      DIST=cost3; 
      MFC=3  ;
    }
  
  _choices_l_l_l_l.createList(input_vertex,reference_vertex);
  _choices_l_l_l_l.putFirst(input_vertex,reference_vertex,MFC);
  
  
  switch(MFC)
    {
    case 1 :  case 2 :
      {
	_choices_l_l_l_l.putLast(input_vertex,reference_vertex,im);
      }
      break;
    case 3 :  case 4 :
      {
	_choices_l_l_l_l.putLast(input_vertex,reference_vertex,jm);
      }
      break;
    case 5 :  case 6 :  case 7 :  case 8 :  case 9 :  case 10 : 
      {
	for (int i=1;i<=T1->getNbChild(input_vertex);i++)
	  {
	    int input_child = T1->child(input_vertex,i);
	    if (input_child==im3)
	      _choices_l_l_l_l.putLast(input_vertex,reference_vertex,jm3);
	    else
	      _choices_l_l_l_l.putLast(input_vertex,reference_vertex,EMPTY_NODE);
	  }
      }
      break;
    default : 	break;
    }
  _d_l_l_l_l.putDBF(input_vertex,reference_vertex,DIST);
  return(DIST);
}


DistanceType MultiscaleMatching::distanceBetweenForest_l_l_v_w(int input_vertex,int reference_vertex) 
{
  int ni=T1->getNbChild(input_vertex);
  int nj=T2->getNbChild(reference_vertex);
  DistanceType cost1,cost2,cost3;
  DistanceType min,DIST;
  int im=0,jm=0,im3=0,jm3=0,MFC=0,mfc=0;
  int i;
  DistanceType dist,dist_l_l_v_w,dist_l_l_v_l,dist_l_l_l_w,dist_v_l,dist_l_w,dist_l_l_l_l;
  
  DIST=MAXDIST;

  min=MAXDIST;
  cost1=getDBF(input_vertex,EMPTY_TREE); 
  for (i=1;i<=ni;i++)
    {
      int input_child=T1->child(input_vertex,i);
      dist_l_l_v_w = _d_l_l_v_w.getDBF(input_child,reference_vertex)-getDBF(input_child,EMPTY_TREE);

      if (dist_l_l_v_w<min)
	{
	  min = dist_l_l_v_w;
	  im  = input_child;
	  mfc = 1;
	}  

      if (!(T1->childIsInComplex(input_vertex,i)))
	{
      dist_l_l_l_w = _d_l_l_l_w.getDBF(input_child,reference_vertex)-getDBF(input_child,EMPTY_TREE);
      dist_v_l = _d_v_l.getDBF(input_child,reference_vertex)-getDBF(input_child,EMPTY_TREE);
	  if (dist_v_l<min)
	    {
	      min = dist_v_l;
	      im  = input_child;
	      mfc = 2;
	    }  
	  if (dist_l_l_l_w<min)
	    {
	      min = dist_l_l_l_w;
	      im  = input_child;
	      mfc = 3;
	    }  
	}
    }
  cost1=cost1+min;
  // On conserve le cout minimum
  if (cost1<DIST)
    {
      DIST=cost1; 
      MFC=mfc; 
    }

  min=MAXDIST;
  cost2=getDBF(EMPTY_TREE,reference_vertex);
  for (i=1;i<=nj;i++)     
    {       
      int reference_child=T2->child(reference_vertex,i);
      dist_l_l_v_w = _d_l_l_v_w.getDBF(input_vertex,reference_child)-getDBF(EMPTY_TREE,reference_child);

      if (dist_l_l_v_w<min)
	{
	  min = dist_l_l_v_w;
	  jm=reference_child;
	  mfc = 4;
	}  


      if (!(T2->childIsInComplex(reference_vertex,i)))
	{
	  dist_l_l_v_l = _d_l_l_v_l.getDBF(input_vertex,reference_child)-getDBF(EMPTY_TREE,reference_child);
	  dist_l_w = _d_l_w.getDBF(input_vertex,reference_child)-getDBF(EMPTY_TREE,reference_child);
	  if (dist_l_w<min) 
	    {
	      min=dist_l_w;
	      jm=reference_child;
	      mfc = 5;
	    }
	  if (dist_l_l_v_l<min)
	    {
	      min = dist_l_l_v_l;
	      jm=reference_child;
	      mfc = 6;
	    }  
	}
    }               
  cost2=cost2+min; 

  if (cost2<DIST) 
    {
      DIST=cost2; 
      MFC=mfc; 
    }
  
  cost3 = MAXDIST;
  NodeList* input_list=new NodeList();
  NodeList* reference_list=new NodeList();
  NEW_MAT(int,choice,ni+1,nj+1)
  IntVector choice_l_l_v_w(ni+1);

  for (int s1=1;s1<=ni;s1++) 
    input_list->push_back(T1->child(input_vertex,s1));

  for (int s2=1;s2<=nj;s2++) 
    reference_list->push_back(T2->child(reference_vertex,s2)); 

 
  _restrMapp_l_l_v_w.make(*input_list,*reference_list);
  _restrMappList_l_l_v_w.resize(ni+nj+3,EMPTY_NODE);
   
  if (ni==0) 
    { 
      _restrMappList_l_l_v_w[1]=2;
      for (i=1;i<=nj;i++)
	{
	  _restrMappList_l_l_v_w[i+1]=1;
	}	  
      cost3=getDBF(EMPTY_TREE,reference_vertex);
   }

  else
    {
      if (nj==0) 
	{ 
	  _restrMappList_l_l_v_w[2]=1;
	  for (i=1;i<=ni;i++)
	    _restrMappList_l_l_v_w[i]=ni+1;	
	  cost3=getDBF(input_vertex,EMPTY_TREE); 
	}
      else
	{
	  DistanceType indelforest =getDBF(input_vertex,EMPTY_TREE)+getDBF(EMPTY_TREE,reference_vertex); 
	  cost3 = indelforest ;
	  for (int s1=1;s1<=ni;s1++) 
	    {
	      int input_child = T1->child(input_vertex,s1);
	      for (int s2=1;s2<=nj;s2++) 
		{
		  int reference_child = T2->child(reference_vertex,s2);
		  dist = MAXDIST;
		  dist_v_l =  _d_v_l.getDBT(input_child,reference_child);
		  dist_l_w =  _d_l_w.getDBT(input_child,reference_child);
		  dist_l_l_l_w =  _d_l_l_l_w.getDBT(input_child,reference_child);
		  dist_l_l_v_w =  _d_l_l_v_w.getDBT(input_child,reference_child);
		  dist_l_l_v_l =  _d_l_l_v_l.getDBT(input_child,reference_child);
		  dist_l_l_l_l =  _d_l_l_l_l.getDBT(input_child,reference_child);
		  if  (T1->childIsInComplex(input_vertex,s1))
		    {
		      if (T2->childIsInComplex(reference_vertex,s2))
			{
			  choice[s1][s2]=1;
			  _restrDistances_l_l_v_w.putDBT(input_child,reference_child,dist_l_l_v_w);
			}
		      else
			{
			  if (dist_l_w<dist)
			    {
			      choice[s1][s2]=2;
			      dist = dist_l_w;
			    }
			  if (dist_l_l_v_l<dist)
			    {
			      choice[s1][s2]=4;
			      dist = dist_l_l_v_l;
			    }
			  if (dist_l_l_l_l<dist)
			    {
			      choice[s1][s2]=6;
			      dist = dist_l_l_l_l;
			    }
			  _restrDistances_l_l_v_w.putDBT(input_child,reference_child,dist);
			}
		    }
		  else
		    {
		      if (T2->childIsInComplex(reference_vertex,s2))
			{
			  if (dist_v_l<dist)
			    {
			      choice[s1][s2]=3;
			      dist = dist_v_l;
			    }
			  if (dist_l_l_l_w<dist)
			    {
			      choice[s1][s2]=5;
			      dist = dist_l_l_l_w;
			    }
			  if (dist_l_l_l_l<dist)
			    {
			      choice[s1][s2]=6;
			      dist = dist_l_l_l_l;
			    }
			  _restrDistances_l_l_v_w.putDBT(input_child,reference_child,dist);
			}
		      else
			{
			  dist =  _d.getDBT(input_child,reference_child);
			  choice[s1][s2]=7;
			  _restrDistances_l_l_v_w.putDBT(input_child,reference_child,dist);
			}
		    }
		}
	    }
	  cost3 = _restrMapp_l_l_v_w.minCostFlow(_restrMappList_l_l_v_w);
	}
    }
		  
		  
  
  if (cost3<=DIST) 
    {
      DIST=cost3; 
      MFC= 7;
    }
  
  _choices_l_l_v_w.createList(input_vertex,reference_vertex);
  _choices_l_l_v_w.putFirst(input_vertex,reference_vertex,MFC);
  
  
  switch(MFC)
    {
    case 1 :   case 2 : case 3 :
      {
	_choices_l_l_v_w.putLast(input_vertex,reference_vertex,im);
      }
      break;
    case 4 :   case 5  : case 6 :
      {
	_choices_l_l_v_w.putLast(input_vertex,reference_vertex,jm);
      }
      break;
    case 7: 
      {
	for (int i=1;i<=ni;i++)
	  {
	    int who = _restrMapp_l_l_v_w.who(_restrMappList_l_l_v_w[i]);
	    if (who==-1)
	      {
		_choices_l_l_v_w.putLast(input_vertex,reference_vertex,-1);
	      }
	    else
	      {
		for (int j=1;j<=nj;j++)
		  {
		    if (T2->child(reference_vertex,j)==who)
		      {
			_choices_l_l_v_w.putLast(input_vertex,reference_vertex,choice[i][j]);
		      }
		  }
	      }
	    _choices_l_l_v_w.putLast(input_vertex,reference_vertex,who);
	  }
      }
      break;
    default : 	break;
    }
  _d_l_l_v_w.putDBF(input_vertex,reference_vertex,DIST);
  DEL_MAT(choice,ni+1)
  return(DIST);
}

DistanceType MultiscaleMatching::distanceBetweenForest_l_l(int input_vertex,int reference_vertex)
{
  
  DistanceType MIN,dist_l_l_v_l,dist_l_l_l_w,dist_l_l_v_w,dist_l_l_l_l;

  dist_l_l_v_l = distanceBetweenForest_l_l_v_l(input_vertex,reference_vertex);
  dist_l_l_l_w = distanceBetweenForest_l_l_l_w(input_vertex,reference_vertex);
  dist_l_l_v_w = distanceBetweenForest_l_l_v_w(input_vertex,reference_vertex);
  dist_l_l_l_l = distanceBetweenForest_l_l_l_l(input_vertex,reference_vertex);
  MIN= D_MIN(D_MIN(D_MIN(dist_l_l_v_l, dist_l_l_l_w), dist_l_l_v_w), dist_l_l_l_l);

  // Il n'est pas nécessaire de retenir les différents choix
  _choices_l_l.createList(input_vertex,reference_vertex);
  _choices_l_l.putFirst(input_vertex,reference_vertex,1);
  _choices_l_l.putFirst(input_vertex,reference_vertex,-1);
  _choices_l_l.putLast(input_vertex,reference_vertex,-1);
  _d_l_l.putDBF(input_vertex,reference_vertex,MIN);
  return(MIN);
}


DistanceType MultiscaleMatching::distanceBetweenForest_v_w(int input_vertex,int reference_vertex)
{
  
  DistanceType MIN=2*MAXDIST,dist_v_w_v_l,dist_v_w_l_w,dist_v_w_v_w,dist_v_w_;
  int MTC = 0,im=-1,jm=-1;

  dist_v_w_v_l = distanceBetweenForest_v_w_v_l(input_vertex,reference_vertex);
  dist_v_w_l_w = distanceBetweenForest_v_w_l_w(input_vertex,reference_vertex);
  dist_v_w_v_w = distanceBetweenForest_v_w_v_w(input_vertex,reference_vertex);
  dist_v_w_ = distanceBetweenForest_v_w_(input_vertex,reference_vertex);

  if (dist_v_w_v_l<MIN)
    {
      MIN = dist_v_w_v_l ;
      MTC = 1 ;
    }
  if (dist_v_w_l_w<MIN)
    {
      MIN = dist_v_w_l_w ;
      MTC = 2;
    }
  if (dist_v_w_v_w<MIN)
    {
      MIN = dist_v_w_v_w ;
      MTC =  3;
    }
  if (dist_v_w_<MIN)
    {
      MIN = dist_v_w_ ;
      MTC =  4;
    }
  _choices_v_w.createList(input_vertex,reference_vertex);
  _choices_v_w.putFirst(input_vertex,reference_vertex,MTC);
  _choices_v_w.putFirst(input_vertex,reference_vertex,-1);
  _choices_v_w.putLast(input_vertex,reference_vertex,-1);
  _d_v_w.putDBF(input_vertex,reference_vertex,MIN);
  return(MIN);

}

DistanceType MultiscaleMatching::distanceBetweenForest(int input_vertex,int reference_vertex)
{
  
  DistanceType MIN=2*MAXDIST,dist_v_l,dist_l_w,dist_v_w,dist_l_l;
  int MTC=0,im=-1, jm=-1;


  dist_v_l = distanceBetweenForest_v_l(input_vertex,reference_vertex);
  dist_l_w = distanceBetweenForest_l_w(input_vertex,reference_vertex);
  dist_v_w = distanceBetweenForest_v_w(input_vertex,reference_vertex);
  dist_l_l = distanceBetweenForest_l_l(input_vertex,reference_vertex);

  if (dist_v_l<MIN)
    {
      MIN = dist_v_l ;
      MTC = 1;
    }
  if (dist_l_w<MIN)
    {
      MIN = dist_l_w ;
      MTC = 2 ;
    }
  if (dist_v_w<MIN)
    {
      MIN = dist_v_w ;
      MTC =  3;
    }
  _choices.createList(input_vertex,reference_vertex);
  _choices.putFirst(input_vertex,reference_vertex,MTC);

  _d.putDBF(input_vertex,reference_vertex,MIN);

  return(MIN);
}

