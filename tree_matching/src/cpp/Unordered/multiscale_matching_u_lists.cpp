#include "multiscale_matching_u.h"



void MultiscaleMatching_U::getList(int input_vertex, int reference_vertex, Sequence* sequence) 
{
  TreeList(input_vertex,reference_vertex,*sequence);
}

// int Lat(ChoiceList* L, int vertex)
// {
//   ChoiceList::iterator begin;
//   begin = L->begin();
//   for (int i=0;i<vertex;i++)
//     begin++;
//   return(*begin);
// }

void MultiscaleMatching_U::ForestList(int input_vertex,int reference_vertex,Sequence& sequence) 
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
    };break;
    default : break;
    }
}

DistanceType MultiscaleMatching_U::match() 
{
  
  DistanceType D=0;
  if ((!T1->isNull())&&(!T2->isNull()))
    {
      for (int input_vertex=T1->getNbVertex()-1;input_vertex>=0;input_vertex--)
	{
	  _d.openDistancesVector(input_vertex);
	  _d_v_l.openDistancesVector(input_vertex);
	  _d_l_w.openDistancesVector(input_vertex);
	  _d_v_w.openDistancesVector(input_vertex);
	  _d_l_l.openDistancesVector(input_vertex);
	  _d_v_w_v_l.openDistancesVector(input_vertex);
	  _d_v_w_l_w.openDistancesVector(input_vertex);
	  _d_v_w_v_w.openDistancesVector(input_vertex);
	  _d_v_w_.openDistancesVector(input_vertex);
	  _d_l_l_v_l.openDistancesVector(input_vertex);
	  _d_l_l_l_w.openDistancesVector(input_vertex);
	  _d_l_l_v_w.openDistancesVector(input_vertex);
	  _d_l_l_l_l.openDistancesVector(input_vertex);
	  _restrDistances_l_l_v_w.openDistancesVector(input_vertex); 
	  _restrDistances_v_w_.openDistancesVector(input_vertex); 
	  _restrDistances_v_w_v_w.openDistancesVector(input_vertex); 
	  _restrDistances.openDistancesVector(input_vertex); 

	  for (int reference_vertex=T2->getNbVertex()-1;reference_vertex>=0;reference_vertex--)
	    {
	      distanceBetweenForest(input_vertex,reference_vertex);
	      distanceBetweenTree(input_vertex,reference_vertex);
	    }
	  
	  for (int i=1;i<=T1->getNbChild(input_vertex);i++)
	    {
	      _d.closeDistancesVector(T1->child(input_vertex,i));
	      _d_v_l.closeDistancesVector(T1->child(input_vertex,i));
	      _d_l_w.closeDistancesVector(T1->child(input_vertex,i));
	      _d_v_w.closeDistancesVector(T1->child(input_vertex,i));
	      _d_l_l.closeDistancesVector(T1->child(input_vertex,i));
	      _d_v_w_v_l.closeDistancesVector(T1->child(input_vertex,i));
	      _d_v_w_l_w.closeDistancesVector(T1->child(input_vertex,i));
	      _d_v_w_v_w.closeDistancesVector(T1->child(input_vertex,i));
	      _d_v_w_.closeDistancesVector(T1->child(input_vertex,i));
	      _d_l_l_v_l.closeDistancesVector(T1->child(input_vertex,i));
	      _d_l_l_l_w.closeDistancesVector(T1->child(input_vertex,i));
	      _d_l_l_v_w.closeDistancesVector(T1->child(input_vertex,i));
	      _d_l_l_l_l.closeDistancesVector(T1->child(input_vertex,i));
	      _restrDistances_l_l_v_w.closeDistancesVector(T1->child(input_vertex,i));
	      _restrDistances_v_w_.closeDistancesVector(T1->child(input_vertex,i)); 
	      _restrDistances_v_w_v_w.closeDistancesVector(T1->child(input_vertex,i)); 
	      _restrDistances.closeDistancesVector(T1->child(input_vertex,i));
	    }
	}

      D=getDBT(0,0);
      cout<<"Distance = "<< D<<endl;
      
    }
  else
    {
      if (T1->isNull())
	{
	  if (!T2->isNull()) {D=_d.referenceTreeFromEmpty(0);};
	}
      else
	{
	  D=_d.inputTreeToEmpty(0);
	}
    }
  return(D);
}



// --------------------------------------------
// Renvoie les distances entre arbres et forets
// --------------------------------------------
DistanceType MultiscaleMatching_U::getDBT(int input_vertex,int reference_vertex) const
{
  return(_d.getDBT(input_vertex,reference_vertex));
}

DistanceType MultiscaleMatching_U::getDBF(int input_vertex,int reference_vertex) const
{
  return(_d.getDBF(input_vertex,reference_vertex));
}





// renvoie le dernier element de la liste de la case node du tableau maintenant les listes d'alignement
int MultiscaleMatching_U::M(int i_node,int r_node)
{
  return(_choices.getList(i_node,r_node)->back());
}		

// renvoie
int MultiscaleMatching_U::M_v_l(int i_node,int r_node) 
{
  return(_choices_v_l.getList(i_node,r_node)->back());
}		

int MultiscaleMatching_U::M_v_w(int i_node,int r_node)
{
  return(_choices_v_w.getList(i_node,r_node)->back());
}		
int MultiscaleMatching_U::M_v_w_v_w(int i_node,int r_node)
{
  return(_choices_v_w_v_w.getList(i_node,r_node)->back());
}		
int MultiscaleMatching_U::M_v_w_v_l(int i_node,int r_node)
{
  return(_choices_v_w_v_l.getList(i_node,r_node)->back());
}		
int MultiscaleMatching_U::M_v_w_l_w(int i_node,int r_node)
{
  return(_choices_v_w_l_w.getList(i_node,r_node)->back());
}		
int MultiscaleMatching_U::M_v_w_(int i_node,int r_node)
{
  return(_choices_v_w_.getList(i_node,r_node)->back());
}		

int MultiscaleMatching_U::M_l_w(int i_node,int r_node)
{
  return(_choices_l_w.getList(i_node,r_node)->back());
}		

int MultiscaleMatching_U::M_l_l(int i_node,int r_node)
{
  return(_choices_l_l.getList(i_node,r_node)->back());
}		
int MultiscaleMatching_U::M_l_l_v_w(int i_node,int r_node)
{
  return(_choices_l_l_v_w.getList(i_node,r_node)->back());
}		
int MultiscaleMatching_U::M_l_l_v_l(int i_node,int r_node)
{
  return(_choices_l_l_v_l.getList(i_node,r_node)->back());
}		
int MultiscaleMatching_U::M_l_l_l_w(int i_node,int r_node)
{
  return(_choices_l_l_l_w.getList(i_node,r_node)->back());
}		
int MultiscaleMatching_U::M_l_l_l_l(int i_node,int r_node)
{
  return(_choices_l_l_l_l.getList(i_node,r_node)->back());
}		

void MultiscaleMatching_U::TreeList(int input_vertex,int reference_vertex,Sequence& sequence)
{
  //cout<<"TreeList "<<input_vertex<<"  "<<reference_vertex<<endl;
  if ((!T1->isNull())&&(!T2->isNull()))
    {
      ChoiceList* L=_choices.getList(input_vertex,reference_vertex);
      int tree_choice=L->front();
      
      switch(tree_choice)
	{
	case 1: 
	  {
	    TreeList_v_l(input_vertex,reference_vertex,sequence);
	  }
	  break;
	case 2: 
	  {
	    TreeList_l_w(input_vertex,reference_vertex,sequence);
	  }
	  break;
	case 3:
	  {
	    TreeList_v_w(input_vertex,reference_vertex,sequence);
	  }
	  break;
	default : break;
	}
    }
}


void MultiscaleMatching_U::TreeList_v_l(int input_vertex,int reference_vertex,Sequence& sequence)
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
	    TreeList_v_l(input_vertex,Lat(L,1),sequence);
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

void MultiscaleMatching_U::TreeList_l_w(int input_vertex,int reference_vertex,Sequence& sequence)
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
	    TreeList_v_w(Lat(L,1),reference_vertex,sequence);
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

void MultiscaleMatching_U::TreeList_v_w_v_l(int input_vertex,int reference_vertex,Sequence& sequence)
{
  //cout<<"TreeList_v_w_v_l"<<input_vertex<<"  "<<reference_vertex<<endl;
  if ((!T1->isNull())&&(!T2->isNull()))
    {
      ChoiceList* L=_choices_v_w_v_l.getList(input_vertex,reference_vertex);
      int tree_choice=L->front();
      
      switch(tree_choice)
	{
	case 1 : 
	  {
	    TreeList_v_w_v_l(Lat(L,1),reference_vertex,sequence);
	  }
	  break;
	case 2 : 
	  {
	    TreeList_v_w_v_l(input_vertex,Lat(L,1),sequence);
	  }
	  break;
	case 3 : 
	  {
	    sequence.append(input_vertex,reference_vertex,_d_v_w_v_l.getCCost(input_vertex,reference_vertex));
	    ForestList_v_w_v_l(input_vertex,reference_vertex,sequence);
	  }
	  break;
	case 4 :
	  {
	    sequence.append(input_vertex,reference_vertex,_d_v_w_v_l.getCCost(input_vertex,reference_vertex));
	    ForestList_l_l_v_l(input_vertex,reference_vertex,sequence);
	  }
	  break;
	default : break;
	}
    }
}

void MultiscaleMatching_U::TreeList_v_w_l_w(int input_vertex,int reference_vertex,Sequence& sequence)
{
  //cout<<"TreeList_v_w_l_w"<<input_vertex<<"  "<<reference_vertex<<endl;
  if ((!T1->isNull())&&(!T2->isNull()))
    {
      ChoiceList* L=_choices_v_w_l_w.getList(input_vertex,reference_vertex);
      int tree_choice=L->front();
      
      switch(tree_choice)
	{
	case 1: 
	  {
	    TreeList_v_w_l_w(Lat(L,1),reference_vertex,sequence);
	  }
	  break;
	case 2: 
	  {
	    TreeList_v_w_l_w(input_vertex,Lat(L,1),sequence);
	  }
	  break;
	case 3 : 
	  {
	    sequence.append(input_vertex,reference_vertex,_d_v_w_l_w.getCCost(input_vertex,reference_vertex));
	    ForestList_v_w_l_w(input_vertex,reference_vertex,sequence);
	  }
	  break;
	case 4 :
	  {
	    sequence.append(input_vertex,reference_vertex,_d_v_w_l_w.getCCost(input_vertex,reference_vertex));
	    ForestList_l_l_l_w(input_vertex,reference_vertex,sequence);
	  }
	  break;
	default : break;
	}
    }
}


void MultiscaleMatching_U::TreeList_v_w_v_w(int input_vertex,int reference_vertex,Sequence& sequence)
{
  //cout<<"TreeList_v_w_v_w"<<input_vertex<<"  "<<reference_vertex<<endl;
  if ((!T1->isNull())&&(!T2->isNull()))
    {
      ChoiceList* L=_choices_v_w_v_w.getList(input_vertex,reference_vertex);
      int tree_choice=L->front();
      
      switch(tree_choice)
	{
	case 1: 
	  {
	    TreeList_v_w_v_w(Lat(L,1),reference_vertex,sequence);
	  }
	  break;
	case 2: 
	  {
	    TreeList_v_w_v_w(input_vertex,Lat(L,1),sequence);
	  }
	  break;
	case 3 : 
	  {
	    sequence.append(input_vertex,reference_vertex,_d_v_w_v_w.getCCost(input_vertex,reference_vertex));
	    ForestList_v_w_v_w(input_vertex,reference_vertex,sequence);
	  }
	  break;
	case 4 :
	  {
	    sequence.append(input_vertex,reference_vertex,_d_v_w_v_w.getCCost(input_vertex,reference_vertex));
	    ForestList_l_l_v_w(input_vertex,reference_vertex,sequence);
	  }
	  break;
	default : break;
	}
    }
}
void MultiscaleMatching_U::TreeList_v_w_(int input_vertex,int reference_vertex,Sequence& sequence)
{
  //cout<<"TreeList_v_w_"<<input_vertex<<"  "<<reference_vertex<<endl;
  if ((!T1->isNull())&&(!T2->isNull()))
    {
      ChoiceList* L=_choices_v_w_.getList(input_vertex,reference_vertex);
      int tree_choice=L->front();
      
      switch(tree_choice)
	{
	case 1: 
	  {
	    TreeList_v_w_(Lat(L,1),reference_vertex,sequence);
	  }
	  break;
	case 2: 
	  {
	    TreeList_v_w_(input_vertex,Lat(L,1),sequence);
	  }
	  break;
	case 3 : 
	  {
	    sequence.append(input_vertex,reference_vertex,_d_v_w_.getCCost(input_vertex,reference_vertex));
	    ForestList_v_w_(input_vertex,reference_vertex,sequence);
	  }
	  break;
	default : break;
	}
    }
}

void MultiscaleMatching_U::TreeList_l_l_v_l(int input_vertex,int reference_vertex,Sequence& sequence)
{
  //cout<<"TreeList_l_l_v_l"<<input_vertex<<"  "<<reference_vertex<<endl;
  if ((!T1->isNull())&&(!T2->isNull()))
    {
      ChoiceList* L=_choices_l_l_v_l.getList(input_vertex,reference_vertex);
      int tree_choice=L->front();
      
      ForestList_l_l_v_l(input_vertex,reference_vertex,sequence);
    }
}

void MultiscaleMatching_U::TreeList_l_l_l_w(int input_vertex,int reference_vertex,Sequence& sequence)
{
  //cout<<"TreeList_l_l_l_w"<<input_vertex<<"  "<<reference_vertex<<endl;
  if ((!T1->isNull())&&(!T2->isNull()))
    {
      ChoiceList* L=_choices_l_l_l_w.getList(input_vertex,reference_vertex);
      int tree_choice=L->front();
      ForestList_l_l_l_w(input_vertex,reference_vertex,sequence);
   }
}


void MultiscaleMatching_U::TreeList_l_l_v_w(int input_vertex,int reference_vertex,Sequence& sequence)
{
  //cout<<"TreeList_l_l_v_w"<<input_vertex<<"  "<<reference_vertex<<endl;
  if ((!T1->isNull())&&(!T2->isNull()))
    {
      ChoiceList* L=_choices_l_l_v_w.getList(input_vertex,reference_vertex);
      int tree_choice=L->front();
      
      ForestList_l_l_v_w(input_vertex,reference_vertex,sequence);
    }
}
void MultiscaleMatching_U::TreeList_l_l_l_l(int input_vertex,int reference_vertex,Sequence& sequence)
{
  //cout<<"TreeList_l_l_l_l"<<input_vertex<<"  "<<reference_vertex<<endl;
  if ((!T1->isNull())&&(!T2->isNull()))
    {
      ChoiceList* L=_choices_l_l_l_l.getList(input_vertex,reference_vertex);
      int tree_choice=L->front();      
      ForestList_l_l_l_l(input_vertex,reference_vertex,sequence);

  }
}
void MultiscaleMatching_U::TreeList_v_w(int input_vertex,int reference_vertex,Sequence& sequence)
{
  //cout<<"TreeList_v_w "<<input_vertex<<"  "<<reference_vertex<<endl;
  if ((!T1->isNull())&&(!T2->isNull()))
    {
      ChoiceList* L=_choices_v_w.getList(input_vertex,reference_vertex);
      int tree_choice=L->front();
      switch(tree_choice)
	{
	case 1 : 
	  {
	    TreeList_v_w_v_l(input_vertex,reference_vertex,sequence);
	  }
	  break;
	case 2 : 
	  {
	    TreeList_v_w_l_w(input_vertex,reference_vertex,sequence);
	  }
	  break;
	case 3:
	  {
	    TreeList_v_w_v_w(input_vertex,reference_vertex,sequence);
	  }
	  break;
	case 4 : 
	  {
	    TreeList_v_w_(input_vertex,reference_vertex,sequence);
	  }
	  break;
	default : break;
	}
    }
}
void MultiscaleMatching_U::TreeList_l_l(int input_vertex,int reference_vertex,Sequence& sequence)
{
  if ((!T1->isNull())&&(!T2->isNull()))
    {
      ChoiceList* L=_choices.getList(input_vertex,reference_vertex);
      int tree_choice=L->front();
      switch(tree_choice)
	{
	case 1 : 
	  {
	    TreeList_l_l_v_l(input_vertex,reference_vertex,sequence);
	  }
	  break;
	case 2 : 
	  {
	    TreeList_l_l_l_w(input_vertex,reference_vertex,sequence);
	  }
	  break;
	case 3:
	  {
	    TreeList_l_l_v_w(input_vertex,reference_vertex,sequence);
	  }
	  break;
	case 4 : 
	  {
	    TreeList_l_l_l_l(input_vertex,reference_vertex,sequence);
	  }
	  break;
	default : break;
	}
    }
}

// Renvoie les matching listes entre forêts

void MultiscaleMatching_U::ForestList_v_l(int input_vertex,int reference_vertex,Sequence& sequence)
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


void MultiscaleMatching_U::ForestList_l_w(int input_vertex,int reference_vertex,Sequence& sequence)
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


void MultiscaleMatching_U::ForestList_v_w_v_l(int input_vertex,int reference_vertex,Sequence& sequence)
{
  //cout<<"ForestList_v_w_v_l"<<input_vertex<<"  "<<reference_vertex<<endl;
  ChoiceList* L=_choices_v_w_v_l.getList(input_vertex,reference_vertex);
  int forest_choice=Lat(L,2);
  
  switch(forest_choice)
    {
    case 1: 
      {
	ForestList_v_w_v_l(Lat(L,3),reference_vertex,sequence);
      }
      break;
    case 2: 
      {
	ForestList_v_w_v_l(input_vertex,Lat(L,3),sequence);
      }
      break;
    case 3: 
      {
	for (int i=1;i<=T1->getNbChild(input_vertex);i++)
	  {
	    int i_node=T1->child(input_vertex,i);
	    int r_node=Lat(L,2+i);
	    if (r_node!=-1) TreeList_v_w_v_l(i_node,r_node,sequence);
	  }
      }
      break;
    default : break;
    }
}


void MultiscaleMatching_U::ForestList_v_w_l_w(int input_vertex,int reference_vertex,Sequence& sequence)
{
  //cout<<"ForestList_v_w_l_w"<<input_vertex<<"  "<<reference_vertex<<endl;
  ChoiceList* L=_choices_v_w_l_w.getList(input_vertex,reference_vertex);
  int forest_choice=Lat(L,2);
  
  switch(forest_choice)
    {
    case 1: 
      {
	ForestList_v_w_l_w(Lat(L,3),reference_vertex,sequence);
      }
      break;
    case 2: 
      {
	ForestList_v_w_l_w(input_vertex,Lat(L,3),sequence);
      }
      break;
    case 3: 
      {
	for (int i=1;i<=T1->getNbChild(input_vertex);i++)
	  {
	    int i_node=T1->child(input_vertex,i);
	    int r_node=Lat(L,2+i);
	    if (r_node!=-1) TreeList_v_w_l_w(i_node,r_node,sequence);
	  }
      }
      break;
    default : break;
    }
}




void MultiscaleMatching_U::ForestList_v_w_(int input_vertex,int reference_vertex,Sequence& sequence)
{
  ChoiceList* L=_choices_v_w_.getList(input_vertex,reference_vertex);
  int forest_choice=Lat(L,2);
  
  switch(forest_choice)
    {
    case 1: 
      {
	ForestList_v_w_(Lat(L,3),reference_vertex,sequence);
      }
      break;
    case 2: 
      {
	ForestList_v_w_(input_vertex,Lat(L,3),sequence);
      }
      break;
    case 3: 
      {
	for (int i=1;i<=T1->getNbChild(input_vertex);i++)
	  {
	    int i_node=T1->child(input_vertex,i);
	    int r_node=Lat(L,2+i);
	    if (r_node!=-1) TreeList_v_w_(i_node,r_node,sequence);
	  }
      }
      break;
    default : break;
    }
}


void MultiscaleMatching_U::ForestList_v_w_v_w(int input_vertex,int reference_vertex,Sequence& sequence)
{
  ChoiceList* L=_choices_v_w_v_w.getList(input_vertex,reference_vertex);
  int forest_choice=Lat(L,2);
  
  switch(forest_choice)
    {
    case 1: 
      {
	ForestList_v_w_v_w(Lat(L,3),reference_vertex,sequence);
      }
      break;
    case 2: 
      {
	ForestList_v_w_v_w(input_vertex,Lat(L,3),sequence);
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
		  if (r_node!=-1) TreeList_v_w_v_w(i_node,r_node,sequence);
		}
		break;
	      case 2: 
		{
		  if (r_node!=-1) TreeList_v_w_(i_node,r_node,sequence);
		}
		break;
	      case 3: 
		{
		  if (r_node!=-1) TreeList_l_l_v_w(i_node,r_node,sequence);
		}
		break;
	      case 4: 
		{
		  if (r_node!=-1) TreeList_l_l_l_l(i_node,r_node,sequence);
		}
		break;
	      case 5: 
		{
		  if (r_node!=-1) TreeList_l_w(i_node,r_node,sequence);
		}
		break;
	      case 6: 
		{
		  if (r_node!=-1) TreeList_v_l(i_node,r_node,sequence);
		}
		break;
	      case 7: 
		{
		  if (r_node!=-1) TreeList_l_l(i_node,r_node,sequence);
		}
		break;
	      case 8: 
		{
		  if (r_node!=-1) TreeList(i_node,r_node,sequence);
		}
		break;
	      default : break;
	      }
	  }
      }
      break;
    default : break;
    }
}

void MultiscaleMatching_U::ForestList_v_w(int input_vertex,int reference_vertex,Sequence& sequence)
{
  ChoiceList* L=_choices_v_w.getList(input_vertex,reference_vertex);
  int forest_choice=Lat(L,2);
  
  switch(forest_choice)
    {
    case 1 : 
      {
	ForestList_v_w_v_l(input_vertex,reference_vertex,sequence);
      }
      break;
    case 2 :
      {
	ForestList_v_w_l_w(input_vertex,reference_vertex,sequence);
      }
      break;
    case 3 :
      {
	ForestList_v_w_v_w(input_vertex,reference_vertex,sequence);
      }
      break;
    case 4 :
      {
	ForestList_v_w_(input_vertex,reference_vertex,sequence);
      }
      break;
    default : break;
    }
}

void MultiscaleMatching_U::ForestList_l_l_v_l(int input_vertex,int reference_vertex,Sequence& sequence)
{
  //cout<<"ForestList_l_l_v_l"<<input_vertex<<"  "<<reference_vertex<<endl;
  ChoiceList* L=_choices_l_l_v_l.getList(input_vertex,reference_vertex);
  int forest_choice=Lat(L,2);
  
  switch(forest_choice)
    {
    case 1: 
      {
	ForestList_l_l_v_l(Lat(L,3),reference_vertex,sequence);
      }
      break;
    case 2: 
      {
	ForestList_l_l_v_l(input_vertex,Lat(L,3),sequence);
      }
      break;
    case 3: 
      {
	ForestList_l_l_v_w(input_vertex,Lat(L,3),sequence);
      }
      break;
    case 4: 
      {
	for (int i=1;i<=T1->getNbChild(input_vertex);i++)
	  {
	    int i_node=T1->child(input_vertex,i);
	    int r_node=Lat(L,2+i);
	    if (r_node!=-1) TreeList_l_l_v_l(i_node,r_node,sequence);
	  }
      }
      break;
    case 5: 
      {
	for (int i=1;i<=T1->getNbChild(input_vertex);i++)
	  {
	    int i_node=T1->child(input_vertex,i);
	    int r_node=Lat(L,2+i);
	    if (r_node!=-1) TreeList_l_l_v_w(i_node,r_node,sequence);
	  }
      }
      break;
    default : break;
    }
}


void MultiscaleMatching_U::ForestList_l_l_l_w(int input_vertex,int reference_vertex,Sequence& sequence)
{
  //cout<<"ForestList_l_l_l_w"<<input_vertex<<"  "<<reference_vertex<<endl;
  ChoiceList* L=_choices_l_l_l_w.getList(input_vertex,reference_vertex);
  int forest_choice=Lat(L,2);
  
  switch(forest_choice)
    {
    case 1: 
      {
	ForestList_l_l_l_w(Lat(L,3),reference_vertex,sequence);
      }
      break;
    case 2: 
      {
	ForestList_l_l_v_w(Lat(L,3),reference_vertex,sequence);
      }
      break;
    case 3: 
      {
	ForestList_l_l_l_w(input_vertex,Lat(L,3),sequence);
      }
      break;
    case 4: 
      {
	for (int i=1;i<=T1->getNbChild(input_vertex);i++)
	  {
	    int i_node=T1->child(input_vertex,i);
	    int r_node=Lat(L,2+i);
	    if (r_node!=-1) TreeList_l_l_l_w(i_node,r_node,sequence);
	  }
      }
      break;
    case 5: 
      {
	for (int i=1;i<=T1->getNbChild(input_vertex);i++)
	  {
	    int i_node=T1->child(input_vertex,i);
	    int r_node=Lat(L,2+i);
	    if (r_node!=-1) TreeList_l_l_v_w(i_node,r_node,sequence);
	  }
      }
      break;
    default : break;
    }
}



void MultiscaleMatching_U::ForestList_l_l_l_l(int input_vertex,int reference_vertex,Sequence& sequence)
{
  //cout<<"ForestList_l_l_l_l"<<input_vertex<<"  "<<reference_vertex<<endl;
  ChoiceList* L=_choices_l_l_l_l.getList(input_vertex,reference_vertex);
  int forest_choice=Lat(L,2);
  
  switch(forest_choice)
    {
    case 1: 
      {
	ForestList_l_l_l_l(Lat(L,3),reference_vertex,sequence);
      }
      break;
    case 2: 
      {
	ForestList_l_l_v_l(Lat(L,3),reference_vertex,sequence);
      }
      break;
    case 3: 
      {
	ForestList_l_l_l_l(input_vertex,Lat(L,3),sequence);
      }
      break;
    case 4: 
      {
	ForestList_l_l_l_w(input_vertex,Lat(L,3),sequence);
      }
      break;
    case 5: 
      {
	for (int i=1;i<=T1->getNbChild(input_vertex);i++)
	  {
	    int i_node=T1->child(input_vertex,i);
	    int r_node=Lat(L,2+i);
	    if (r_node!=-1) TreeList_l_l_l_l(i_node,r_node,sequence);
	  }
      }
      break;
    case 6: 
      {
	for (int i=1;i<=T1->getNbChild(input_vertex);i++)
	  {
	    int i_node=T1->child(input_vertex,i);
	    int r_node=Lat(L,2+i);
	    if (r_node!=-1) TreeList_l_w(i_node,r_node,sequence);
	  }
      }
      break;
    case 7: 
      {
	for (int i=1;i<=T1->getNbChild(input_vertex);i++)
	  {
	    int i_node=T1->child(input_vertex,i);
	    int r_node=Lat(L,2+i);
	    if (r_node!=-1) TreeList_v_l(i_node,r_node,sequence);
	  }
      }
      break;
    case 8: 
      {
	for (int i=1;i<=T1->getNbChild(input_vertex);i++)
	  {
	    int i_node=T1->child(input_vertex,i);
	    int r_node=Lat(L,2+i);
	    if (r_node!=-1) TreeList_l_l_v_l(i_node,r_node,sequence);
	  }
      }
      break;
    case 9: 
      {
	for (int i=1;i<=T1->getNbChild(input_vertex);i++)
	  {
	    int i_node=T1->child(input_vertex,i);
	    int r_node=Lat(L,2+i);
	    if (r_node!=-1) TreeList_l_l_l_w(i_node,r_node,sequence);
	  }
      }
      break;
    case 10: 
      {
	for (int i=1;i<=T1->getNbChild(input_vertex);i++)
	  {
	    int i_node=T1->child(input_vertex,i);
	    int r_node=Lat(L,2+i);
	    if (r_node!=-1) TreeList_l_l_l_w(i_node,r_node,sequence);
	  }
      }
      break;
    default : break;
    }
}



void MultiscaleMatching_U::ForestList_l_l_v_w(int input_vertex,int reference_vertex,Sequence& sequence)
{
  ChoiceList* L=_choices_l_l_v_w.getList(input_vertex,reference_vertex);
  int forest_choice=Lat(L,2);
  
  
  switch(forest_choice)
    {
    case 1: 
      {
	ForestList_l_l_v_w(Lat(L,3),reference_vertex,sequence);
      }
      break;
    case 2: 
      {
	ForestList_v_l(Lat(L,3),reference_vertex,sequence);
      }
      break;
    case 3: 
      {
	ForestList_l_l_l_w(Lat(L,3),reference_vertex,sequence);
      }
      break;
    case 4: 
      {
	ForestList_l_l_v_w(input_vertex,Lat(L,3),sequence);
      }
      break;
    case 5: 
      {
	ForestList_l_w(input_vertex,Lat(L,3),sequence);
      }
      break;
    case 6: 
      {
	ForestList_l_l_v_l(input_vertex,Lat(L,3),sequence);
      }
      break;
    case 7: 
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
		  if (r_node!=-1) TreeList_l_l_v_w(i_node,r_node,sequence);
		}
		break;
	      case 2: 
		{
		  if (r_node!=-1) TreeList_l_w(i_node,r_node,sequence);
		}
		break;
	      case 3: 
		{
		  if (r_node!=-1) TreeList_v_l(i_node,r_node,sequence);
		}
		break;
	      case 4: 
		{
		  if (r_node!=-1) TreeList_l_l_v_l(i_node,r_node,sequence);
		}
		break;
	      case 5: 
		{
		  if (r_node!=-1) TreeList_l_l_l_w(i_node,r_node,sequence);
		}
		break;
	      case 6: 
		{
		  if (r_node!=-1) TreeList_l_l_l_l(i_node,r_node,sequence);
		}
		break;
	      case 7: 
		{
		  if (r_node!=-1) TreeList(i_node,r_node,sequence);
		}
		break;
	      default : break;
	      }
	  }
      }
      break;
    default : break;
    }
}




