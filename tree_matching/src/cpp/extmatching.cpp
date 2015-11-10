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


#include "extmatching.h"

const int EMPTY_FOREST = -1;

ExtendedMatching::ExtendedMatching(TreeGraph& input_tree,TreeGraph& reference_tree,NodeCost& node_distance)
{
  T1=&input_tree;
  T2=&reference_tree;
  ND=&node_distance;
  _distances.make(*T1,*T2,node_distance);
  _choices.resize(T1->getNbVertex(),T2->getNbVertex());
  _restrMappIMFAndRMF.link(I_MAX(T1->getDegree(),T2->getDegree()),*_distances.getDistanceTable());
  _restrMappIPFAndRPF.link(I_MAX(T1->getDegree(),T2->getDegree()),*_distances.getDistanceTable());
  _restrMappIMFAndRF.link(I_MAX(T1->getDegree(),T2->getDegree()),*_distances.getDistanceTable());
  _restrMappIFAndRMF.link(I_MAX(T1->getDegree(),T2->getDegree()),*_distances.getDistanceTable());
}

ExtendedMatching::~ExtendedMatching()
{
}

DistanceType ExtendedMatching::distanceBetweenTree(int input_vertex,int reference_vertex)
{
  int ni=T1->getNbChild(input_vertex);
  int nj=T2->getNbChild(reference_vertex);
  DistanceType cost1,cost2,cost3,dist1,dist2;
  DistanceType min,MIN=2*MAXDIST;
  int im=0,jm=0,MTC=0;
  int i;

//--------------------------------------------------------------------------------------------------------------
//Case 1 : We search the reference_tree as a subtree of the input_tree
//--------------------------------------------------------------------------------------------------------------

  min=MAXDIST;
  cost1=getDBT(input_vertex,EMPTY_TREE);
  for (i=1;i<=ni;i++)
  {
    int input_child=T1->child(input_vertex,i);
    dist1=getDBT(input_child,reference_vertex)-getDBT(input_child,EMPTY_TREE);
    if (dist1<min) { min=dist1; im=input_child; };
  }
  cost1=cost1+min;

  if (cost1<MIN) { MIN=cost1; MTC=1; };

//--------------------------------------------------------------------------------------------------------------
//Case2 : We search the input_tree as a subtree of the reference_tree
//--------------------------------------------------------------------------------------------------------------

  min=MAXDIST;
  cost2=getDBT(EMPTY_TREE,reference_vertex);
  for (i=1;i<=nj;i++)
  {
    int reference_child=T2->child(reference_vertex,i);
    dist2=getDBT(input_vertex,reference_child)-getDBT(EMPTY_TREE,reference_child);
    if (dist2<min) { min=dist2;jm=reference_child; };
  }
  cost2=cost2+min;

  if (cost2<MIN) { MIN=cost2; MTC=2; };

//--------------------------------------------------------------------------------------------------------------
//Case3 : We evaluate the matching between the input_forest and the reference_forest
//--------------------------------------------------------------------------------------------------------------

  cost3=getDBF(input_vertex,reference_vertex);
  cost3=cost3+_distances.getCCost(input_vertex,reference_vertex);
  if (cost3<MIN) { MIN=cost3; MTC=3;  };

//--------------------------------------------------------------------------------------------------------------
// We maintain the matching lists
//--------------------------------------------------------------------------------------------------------------
  switch (MTC)
  {
  case 1 :
    {
      _choices.putFirst(input_vertex,reference_vertex,im);
      _choices.putLast(input_vertex,reference_vertex,-1);
    };break;
  case 2 :
    {
      _choices.putFirst(input_vertex,reference_vertex,jm);
      _choices.putLast(input_vertex,reference_vertex,M(input_vertex,jm));
    };break;
  case 3 :
    {
      _choices.putFirst(input_vertex,reference_vertex,-1);
      _choices.putLast(input_vertex,reference_vertex,reference_vertex);
    };break;
    default :   assert(0);break;
  }
  _choices.putFirst(input_vertex,reference_vertex,MTC);

  _distances.putDBT(input_vertex,reference_vertex,MIN);
  return(MIN);
}




DistanceType ExtendedMatching::distanceBetweenForest(int input_vertex,int reference_vertex)
{
  int if_nb_tree = T1->getNbChild(input_vertex);
  int rf_nb_tree = T2->getNbChild(reference_vertex);
  DistanceType dist1,dist2,dist3;
  int ifm,jfm,ifm_mf,jfm_mf,ifm_pf,jfm_pf,ifm_ifrmf,jfm_ifrmf,ifm_imfrf,jfm_imfrf;
  DistanceType tmp_dist1,tmp_dist2,tmp_db_pf,tmp_db_mf,tmp_db_imfrf,tmp_db_ifrmf;
  DistanceType db_pf,db_mf,db_imfrf,db_ifrmf;
  int mc_f,mc_pf,mc_mf,mc_imfrf,mc_ifrmf;
  DistanceType db_f = MAXDIST;
  int son;

//--------------------------------------------------------------------------------------------------------------
//Case 1 : We search the reference_forest as a sub-forest of the input_forest
//--------------------------------------------------------------------------------------------------------------
  dist1    = MAXDIST;
  db_pf    = MAXDIST;
  db_mf    = MAXDIST;
  db_imfrf = MAXDIST;
  db_ifrmf = MAXDIST;

  for (son=1;son<=if_nb_tree;son++)
  {
    int input_child=T1->child(input_vertex,son);
    if (sameComplex(T1->getNode(input_child),T1->getNode(input_vertex)))
    {

      tmp_dist1 = getDBPF(input_child,reference_vertex)-getDBPF(input_child,EMPTY_FOREST);
      tmp_dist1 = tmp_dist1+getDBMF(input_child,reference_vertex)-getDBMF(input_child,EMPTY_FOREST);
      if (tmp_dist1<dist1)
      {
        dist1=tmp_dist1;
        ifm=input_child;
      };

      tmp_db_pf = getDBPF(input_child,reference_vertex)-getDBPF(input_child,EMPTY_FOREST)+getDBPF(input_vertex,EMPTY_FOREST);
      if (tmp_db_pf<db_pf)
      {
        db_pf=tmp_db_pf;
        mc_pf=1;
        ifm_pf=input_child;
      };

      tmp_db_ifrmf=getDBF(input_vertex,EMPTY_FOREST)-getDBMF(input_child,EMPTY_FOREST)+getDBMF(input_child,reference_vertex);
      if (tmp_db_ifrmf<db_ifrmf)
      {
        db_ifrmf=tmp_db_ifrmf;
        mc_ifrmf=1;
        ifm_ifrmf=input_child;
      }
    }
    else
    {
      tmp_dist1 = getDBIFRMF(input_child,reference_vertex)-getDBMF(input_child,EMPTY_FOREST);
      tmp_dist1 = tmp_dist1+getDBPF(EMPTY_FOREST,reference_vertex);
      if (tmp_dist1<dist1)
      {
        dist1=tmp_dist1;
        ifm=input_child;
      };

      tmp_db_mf = getDBIFRMF(input_child,reference_vertex)-getDBF(input_child,EMPTY_FOREST)+getDBMF(input_vertex,EMPTY_FOREST);
      if (tmp_db_mf<db_mf)
      {
        db_mf=tmp_db_mf;
        mc_mf=1;
        ifm_mf=input_child;
      };

      tmp_db_ifrmf=getDBF(input_child,EMPTY_FOREST)-getDBF(input_child,EMPTY_FOREST)+getDBIFRMF(input_child,reference_vertex);
      if (tmp_db_ifrmf<db_ifrmf)
      {
        db_ifrmf=tmp_db_ifrmf;
        mc_ifrmf=1;
        ifm_ifrmf=input_child;
      };

      tmp_db_imfrf=getDBMF(input_child,EMPTY_FOREST)+getDBPF(EMPTY_FOREST,reference_vertex)+getDBMF(input_child,reference_vertex)-getDBMF(input_child,EMPTY_FOREST);
      if (tmp_db_imfrf<db_imfrf)
      {
        db_imfrf=tmp_db_imfrf;
        mc_imfrf=1;
        ifm_imfrf=input_child;
      };

    }
  }

  dist1=getDBF(input_vertex,EMPTY_FOREST)+dist1;

  if (dist1<=db_f)
  {
    db_f = dist1;
    mc_f = 1;
  };

//--------------------------------------------------------------------------------------------------------------
//Case 2 : We search the input_forest as a sub-forest of the reference_forest
//--------------------------------------------------------------------------------------------------------------
  dist2 = MAXDIST;

  for (son=1;son<=rf_nb_tree;son++)
  {
    int reference_child=T2->child(reference_vertex,son);
    if (sameComplex(T2->getNode(reference_child),T2->getNode(reference_vertex)))
    {
      tmp_dist2 = getDBPF(input_vertex,reference_child)-getDBPF(EMPTY_FOREST,reference_child);
      tmp_dist2 = tmp_dist2+getDBMF(input_vertex,reference_child)-getDBMF(EMPTY_FOREST,reference_child);
      if (tmp_dist2<dist2)
      {
        dist2=tmp_dist2;
        jfm=reference_child;
      }

      tmp_db_pf = getDBPF(input_vertex,reference_child)-getDBPF(EMPTY_FOREST,reference_child)+getDBPF(EMPTY_FOREST,reference_vertex);
      if (tmp_db_pf<db_pf)
      {
        db_pf=tmp_db_pf;
        mc_pf=2;
        jfm_pf=reference_child;
      }

      tmp_db_imfrf=getDBF(EMPTY_FOREST,reference_vertex)-getDBMF(EMPTY_FOREST,reference_child)+getDBMF(input_vertex,reference_child);
      if (tmp_db_imfrf<db_imfrf)
      {
        db_imfrf=tmp_db_imfrf;
        mc_imfrf=2;
        jfm_imfrf=reference_child;
      };

    }
    else
    {
      tmp_dist2 = getDBIMFRF(input_vertex,reference_child)-getDBF(EMPTY_FOREST,reference_child);
      tmp_dist2 = tmp_dist2 + getDBPF(input_vertex,EMPTY_FOREST);
      if (tmp_dist2<dist2)
      {
        dist2 = tmp_dist2;
        jfm=reference_child;
      };

      tmp_db_mf = getDBIFRMF(input_vertex,reference_child)-getDBF(EMPTY_FOREST,reference_child)+getDBMF(EMPTY_FOREST,reference_vertex);
      if (tmp_db_mf<db_mf)
      {
        db_mf=tmp_db_mf;
        mc_mf=2;
        jfm_mf=reference_child;
      };

      tmp_db_imfrf = getDBF(EMPTY_FOREST,reference_vertex)-getDBF(EMPTY_FOREST,reference_child)+getDBIMFRF(input_vertex,reference_child);
      if (tmp_db_imfrf<db_imfrf)
      {
        db_imfrf=tmp_db_imfrf;
        mc_imfrf=2;
        jfm_imfrf=reference_child;
      };

      tmp_db_ifrmf = getDBPF(input_vertex,EMPTY_FOREST)+getDBPF(EMPTY_FOREST,reference_vertex)+getDBMF(input_vertex,reference_child)-getDBMF(EMPTY_FOREST,reference_child);
      if (tmp_db_ifrmf<db_ifrmf)
      {
        db_ifrmf=tmp_db_ifrmf;
        mc_ifrmf=2;
        jfm_ifrmf=reference_child;
      };
    }
  }

  dist2 = getDBF(EMPTY_FOREST,reference_vertex) + dist2;

  if (dist2<=db_f)
  {
    db_f = dist2;
    mc_f = 2;
  };

//--------------------------------------------------------------------------------------------------------------
//Case 3 : We evaluate the restricted mapping between the input_forest and the reference_forest
//--------------------------------------------------------------------------------------------------------------
  NodeList* ipf_list = new NodeList();
  NodeList* imf_list = new NodeList();
  NodeList* rpf_list = new NodeList();
  NodeList* rmf_list = new NodeList();
  NodeList* if_list  = new NodeList();
  NodeList* rf_list  = new NodeList();

  for (son=1;son<=if_nb_tree;son++)
  {
    int input_son=T1->child(input_vertex,son);
    if_list->push_back(input_son);
    if (sameComplex(T1->getNode(input_son),T1->getNode(input_vertex)))
    {
      ipf_list->push_back(input_son);
    }
    else
    {
      imf_list->push_back(input_son);
    }
  }

  for (son=1;son<=rf_nb_tree;son++)
  {
    int reference_son=T2->child(reference_vertex,son);
    rf_list->push_back(reference_son);
    if (sameComplex(T2->getNode(reference_son),T2->getNode(reference_vertex)))
    {
      rpf_list->push_back(reference_son);
    }
    else
    {
      rmf_list->push_back(reference_son);
    }
  }

  int ipf_nb_tree = ipf_list->size();
  int imf_nb_tree = imf_list->size();
  int rpf_nb_tree = rpf_list->size();
  int rmf_nb_tree = rmf_list->size();

  //--------------------------------------------------------------------------------------
  // MATCHING BETWEEN THE INPUT MINUS FOREST AND THE  REFERENCE FOREST
  //--------------------------------------------------------------------------------------

  _restrMappIMFAndRF.make(*imf_list,*rf_list);
  _restrMappIMFAndRFList.resize(imf_nb_tree+rf_nb_tree+3, EMPTY_NODE);


  int i;
  // THE INPUT MINUS FOREST IS EMPTY
  // All the reference vertices are paired with empty
  if (imf_nb_tree==0)
  {
     _restrMappIMFAndRFList[1]=2;
    for (i=1;i<=rf_nb_tree;i++) { _restrMappIMFAndRFList[i+1]=1; };
    tmp_db_imfrf=getDBF(EMPTY_FOREST,reference_vertex);
  }
  else
  {
    // THE REFERENCE FOREST IS EMPTY
    // All the input vertices are paired with empty
    if (rf_nb_tree==0)
    {
      _restrMappIMFAndRFList[2]=1;
      for (i=1;i<=imf_nb_tree;i++) { _restrMappIMFAndRFList[i]=imf_nb_tree+1; };
      tmp_db_imfrf=getDBF(input_vertex,EMPTY_FOREST);
    }
    else
    {
      // BOTH FOREST ARE NOT EMPTY
      // A retricted mapping must be calculated
      tmp_db_imfrf = _restrMappIMFAndRF.minCostFlow(_restrMappIMFAndRFList);
    }
  }

  if (tmp_db_imfrf<db_imfrf)
  {
    db_imfrf=tmp_db_imfrf;
    mc_imfrf=3;
  };

  _distances.putDBIMFRF(input_vertex,reference_vertex,db_imfrf);

  //--------------------------------------------------------------------------------------
  // MATCHING BETWEEN THE INPUT FOREST AND THE REFERENCE MINUS FOREST
  //--------------------------------------------------------------------------------------

  _restrMappIFAndRMF.make(*if_list,*rmf_list);
  _restrMappIFAndRMFList.resize(if_nb_tree+rmf_nb_tree+3, EMPTY_NODE);


  // THE INPUT FOREST IS EMPTY
  // All the reference vertices are paired with empty
  if (if_nb_tree==0)
  {
    _restrMappIFAndRMFList[1]=2;
    for (i=1;i<=rmf_nb_tree;i++) { _restrMappIFAndRMFList[i+1]=1; };
    tmp_db_ifrmf=getDBMF(EMPTY_FOREST,reference_vertex);
  }
  else
  {
    // THE REFERENCE MINUS FOREST IS EMPTY
    // All the input vertices are paired with empty
    if (rmf_nb_tree==0)
    {
      _restrMappIFAndRMFList[2]=1;
      for (i=1;i<=if_nb_tree;i++) { _restrMappIFAndRMFList[i]=if_nb_tree+1; };
      tmp_db_ifrmf=getDBMF(input_vertex,EMPTY_FOREST);
    }
    else
    {
      // BOTH FOREST ARE NOT EMPTY
      // A retricted mapping must be calculated
      tmp_db_ifrmf=_restrMappIFAndRMF.minCostFlow(_restrMappIFAndRMFList);
    }
  }

  if (tmp_db_ifrmf<db_ifrmf)
  {
    db_ifrmf=tmp_db_ifrmf;
    mc_ifrmf=3;
  };

  _distances.putDBIFRMF(input_vertex,reference_vertex,db_ifrmf);

  //--------------------------------------------------------------------------------------
  // MATCHING BETWEEN THE INPUT PLUS FOREST AND THE REFERENCE PLUS FOREST
  //--------------------------------------------------------------------------------------

  _restrMappIPFAndRPF.make(*ipf_list,*rpf_list);
  _restrMappIPFAndRPFList.resize(ipf_nb_tree+rpf_nb_tree+3, EMPTY_NODE);


  if (ipf_nb_tree==0)
  {
    // THE INPUT PLUS FOREST IS EMPTY
    // All the reference vertices are paired with empty
    _restrMappIPFAndRPFList[1]=2;
    for (i=1;i<=rpf_nb_tree;i++) { _restrMappIPFAndRPFList[i+1]=1; };
    tmp_db_pf = getDBPF(EMPTY_FOREST,reference_vertex);
   }
  else
  {
    // THE REFERENCE PLUS FOREST IS EMPTY
    // All the input vertices are paired with empty
    if (rpf_nb_tree==0)
    {
      _restrMappIPFAndRPFList[2]=1;
      for (i=1;i<=ipf_nb_tree;i++) { _restrMappIPFAndRPFList[i]=ipf_nb_tree+1; };
      tmp_db_pf = getDBPF(input_vertex,EMPTY_FOREST);
    }
    else
    {
      //BOTH FOREST ARE NOT EMPTY
      // A retricted mapping must be calculated
      tmp_db_pf = _restrMappIPFAndRPF.minCostFlow(_restrMappIPFAndRPFList);
    }
  }

  if (tmp_db_pf<db_pf)
  {
    db_pf=tmp_db_pf;
    mc_pf=3;
  };

  _distances.putDBPF(input_vertex,reference_vertex,db_pf);

  //--------------------------------------------------------------------------------------
  // MATCHING BETWEEN THE INPUT MINUS FOREST AND THE  REFERENCE FOREST
  //--------------------------------------------------------------------------------------

  _restrMappIMFAndRMF.make(*imf_list,*rmf_list);
  _restrMappIMFAndRMFList.resize(imf_nb_tree + rmf_nb_tree + 3, EMPTY_NODE);


  // THE INPUT MINUS FOREST IS EMPTY
  // All the reference vertices are paired with empty
  if (imf_nb_tree==0)
  {
    _restrMappIMFAndRMFList[1]=2;
    for (i=1;i<=rmf_nb_tree;i++) { _restrMappIMFAndRMFList[i+1]=1; };
    tmp_db_mf=getDBMF(EMPTY_FOREST,reference_vertex);
  }
  else
  {
    // THE REFERENCE MINUS FOREST IS EMPTY
    // All the input vertices are paired with empty
    if (rmf_nb_tree==0)
    {
      _restrMappIMFAndRMFList[2]=1;
      for (i=1;i<=imf_nb_tree;i++) { _restrMappIMFAndRMFList[i]=imf_nb_tree+1; };
      tmp_db_mf=getDBMF(input_vertex,EMPTY_FOREST);
    }
    else
    {
      // BOTH FOREST ARE NOT EMPTY
      // A retricted mapping must be calculated
      tmp_db_mf=_restrMappIMFAndRMF.minCostFlow(_restrMappIMFAndRMFList);
    }
  }

  if (tmp_db_mf<db_mf)
  {
    db_mf=tmp_db_mf;
    mc_mf=3;
  };

  _distances.putDBMF(input_vertex,reference_vertex,db_mf);

  dist3 = db_pf + db_mf;

  if (dist3<=db_f)
  {
    db_f = dist3;
    mc_f = 3;
  };

//--------------------------------------------------------------------------------------------------------------
//We maintain the matching lists
//--------------------------------------------------------------------------------------------------------------
  _choices.createList(input_vertex,reference_vertex);
  _choices.putIPFRPFLast(input_vertex,reference_vertex,mc_pf);
  switch(mc_pf)
  {
  case 1 :
    {
     _choices.putIPFRPFLast(input_vertex,reference_vertex,ifm_pf);
    };
    break;
  case 2 :
    {
      _choices.putIPFRPFLast(input_vertex,reference_vertex,jfm_pf);
    };
    break;
  case 3 :
    {
     int s;
      for (s=1;s<=ipf_list->size();s++)
      {
        _choices.putIPFRPFLast(input_vertex,reference_vertex,_restrMappIPFAndRPF.who(_restrMappIPFAndRPFList[s]));
      };
    };
    break;
    default : {assert(0);};  break;
  };


 _choices.putIMFRMFLast(input_vertex,reference_vertex,mc_mf);
   switch(mc_mf)
  {
  case 1 :
    {
      _choices.putIMFRMFLast(input_vertex,reference_vertex,ifm_mf);
    };
    break;
  case 2 :
    {
      _choices.putIMFRMFLast(input_vertex,reference_vertex,jfm_mf);
    };
    break;
  case 3 :
    {
      int s;
      for (s=1;s<=imf_list->size();s++)
      {
        _choices.putIMFRMFLast(input_vertex,reference_vertex,_restrMappIMFAndRMF.who(_restrMappIMFAndRMFList[s]));
      };
    };
    break;
    default : {assert(0);}; break;
  };

  _choices.putIFRMFLast(input_vertex,reference_vertex,mc_ifrmf);
  switch(mc_ifrmf)
  {
  case 1 :
    {
      _choices.putIFRMFLast(input_vertex,reference_vertex,ifm_ifrmf);
    };
    break;
  case 2 :
    {
      _choices.putIFRMFLast(input_vertex,reference_vertex,jfm_ifrmf);
    };
    break;
  case 3 :
    {
      int s;
      for (s=1;s<=if_list->size();s++)
      {
        _choices.putIFRMFLast(input_vertex,reference_vertex,_restrMappIFAndRMF.who(_restrMappIFAndRMFList[s]));
      };
    };
    break;
    default :{assert(0);}; break;
  };

  _choices.putIMFRFLast(input_vertex,reference_vertex,mc_imfrf);
   switch(mc_imfrf)
  {
  case 1 :
    {
      _choices.putIMFRFLast(input_vertex,reference_vertex,ifm_imfrf);
    };
    break;
  case 2 :
    {
      _choices.putIMFRFLast(input_vertex,reference_vertex,jfm_imfrf);
    };
    break;
  case 3 :
    {
      int s;
      for (s=1;s<=imf_list->size();s++)
      {
        _choices.putIMFRFLast(input_vertex,reference_vertex,_restrMappIMFAndRF.who(_restrMappIMFAndRFList[s]));
      };
    };
    break;
    default :{ assert(0);}; break;
  };

  _choices.putFirst(input_vertex,reference_vertex,mc_f);
  switch(mc_f)
  {
  case 1 :
    {
      _choices.putLast(input_vertex,reference_vertex,ifm);
   };
    break;
  case 2 :
    {
      _choices.putLast(input_vertex,reference_vertex,jfm);
    };
    break;
  case 3 :
    {
      int s;
      for (s=1;s<=ipf_list->size();s++)
      {
        _choices.putLast(input_vertex,reference_vertex,_restrMappIPFAndRPF.who(_restrMappIPFAndRPFList[s]));
      };

      for (s=1;s<imf_list->size();s++)
      {
        _choices.putLast(input_vertex,reference_vertex,_restrMappIMFAndRMF.who(_restrMappIMFAndRMFList[s]));
      };

    };
    break;
    default : {assert(0);}; break;
  };

  delete (NodeList*) if_list;
  delete (NodeList*) rf_list;
  delete (NodeList*) ipf_list;
  delete (NodeList*) rpf_list;
  delete (NodeList*) imf_list;
  delete (NodeList*) rmf_list;

  _distances.putDBF(input_vertex,reference_vertex,db_f);
  return(db_f);

}

void ExtendedMatching::getList(int input_vertex, int reference_vertex, Sequence* sequence)
{
  TreeList(input_vertex,reference_vertex,*sequence);
}

int ExtendedMatching::Lat(ChoiceList* L, int vertex)
{
  ChoiceList::iterator begin;
  begin = L->begin();
  for (int i=0;i<vertex;i++)
    begin++;
  return(*begin);
}



void ExtendedMatching::TreeList(int input_vertex,int reference_vertex,Sequence& sequence)
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
      };break;
    case 2:
      {
        TreeList(input_vertex,Lat(L,1),sequence);
      };break;
    case 3:
      {
        sequence.append(input_vertex,reference_vertex,_distances.getCCost(input_vertex,reference_vertex));
        ForestList(input_vertex,reference_vertex,sequence);
      };break;
      default : break;
    }
  }
}

void ExtendedMatching::ForestList(int input_vertex,int reference_vertex,Sequence& sequence)
{
  ChoiceList* L=_choices.getList(input_vertex,reference_vertex);
  int forest_choice=Lat(L,2);
  switch(forest_choice)
  {
  case 1:
    {
      int input_child=Lat(L,3);
      if (sameComplex(T1->getNode(input_child),T1->getNode(input_vertex)))
      {
        pfList(input_child,reference_vertex,sequence);
        mfList(input_child,reference_vertex,sequence);
      }
      else
      {
        imfrfList(input_child,reference_vertex,sequence);
      };
    };break;
  case 2:
    {
      int reference_child=Lat(L,3);
      if (sameComplex(T2->getNode(reference_child),T2->getNode(reference_vertex)))
      {
        pfList(input_vertex,reference_child,sequence);
        mfList(input_vertex,reference_child,sequence);
      }
      else
      {
        ifrmfList(input_vertex,reference_child,sequence);
      }
    };break;
  case 3:
    {
      int i;
      int child_index=1;
      for (i=1;i<=T1->getNbChild(input_vertex);i++)
      {
        int input_child = T1->child(input_vertex,i);
        if (sameComplex(T1->getNode(input_child),T1->getNode(input_vertex)))
        {
          int reference_child = Lat(L,2+child_index);
          child_index++;
          if (reference_child != EMPTY_NODE)
          {
            TreeList(input_child,reference_child,sequence);
          }
        }
      };

      for (i=1;i<=T1->getNbChild(input_vertex);i++)
      {
        int input_child = T1->child(input_vertex,i);
        if (!sameComplex(T1->getNode(input_child),T1->getNode(input_vertex)))
        {
          int reference_child = Lat(L,2+child_index);
          child_index++;
          if (reference_child != EMPTY_NODE)
          {
            TreeList(input_child,reference_child,sequence);
          }
        }
      };

    };break;
    default : break;
  }
}

void ExtendedMatching::pfList(int input_vertex,int reference_vertex,Sequence& sequence)
{

  ChoiceList* L=_choices.getIPFRPFList(input_vertex,reference_vertex);
  int forest_choice=Lat(L,0);
  switch(forest_choice)
  {
  case 1:
    {
      pfList(Lat(L,1),reference_vertex,sequence);
    };break;
  case 2:
    {
      pfList(input_vertex,Lat(L,1),sequence);
    };break;
  case 3:
    {
      for (int i=1;i<=T1->getNbChild(input_vertex);i++)
      {
        int child_index=1;
        int input_child = T1->child(input_vertex,i);
        if (sameComplex(T1->getNode(input_child),T1->getNode(input_vertex)))
        {
          int reference_child=Lat(L,child_index++);
          if (reference_child != EMPTY_NODE)
          {
            TreeList(input_child,reference_child,sequence);
          }
        }
      }
    };break;
    default : break;
  }
}

void ExtendedMatching::mfList(int input_vertex,int reference_vertex,Sequence& sequence)
{
  ChoiceList* L=_choices.getIMFRMFList(input_vertex,reference_vertex);
  int forest_choice=Lat(L,0);
  switch(forest_choice)
  {
  case 1:
    {
      mfList(Lat(L,1),reference_vertex,sequence);
    };break;
  case 2:
    {
      mfList(input_vertex,Lat(L,1),sequence);
    };break;
  case 3:
    {
      for (int i=1;i<=T1->getNbChild(input_vertex);i++)
      {
        int child_index=1;
        int input_child = T1->child(input_vertex,i);
        if (!sameComplex(T1->getNode(input_child),T1->getNode(input_vertex)))
        {
          int reference_child=Lat(L,child_index++);
          if (reference_child != EMPTY_NODE)
          {
            TreeList(input_child,reference_child,sequence);
          }
        }
      }
    };break;
    default : break;
  }
};

void ExtendedMatching::ifrmfList(int input_vertex,int reference_vertex,Sequence& sequence)
{
  ChoiceList* L=_choices.getIFRMFList(input_vertex,reference_vertex);
  int forest_choice=Lat(L,0);
  switch(forest_choice)
  {
  case 1:
    {
      ifrmfList(Lat(L,1),reference_vertex,sequence);
    };break;
  case 2:
    {
      ifrmfList(input_vertex,Lat(L,1),sequence);
    };break;
  case 3:
    {
      for (int i=1;i<=T1->getNbChild(input_vertex);i++)
      {
        int input_child = T1->child(input_vertex,i);
        int reference_child=Lat(L,i);
        if (reference_child != EMPTY_NODE)
        {
          TreeList(input_child,reference_child,sequence);
        }
      }
    };break;
    default : break;
  }
}

void ExtendedMatching::imfrfList(int input_vertex,int reference_vertex,Sequence& sequence)
{
  ChoiceList* L=_choices.getIMFRFList(input_vertex,reference_vertex);
  int forest_choice=Lat(L,0);
  switch(forest_choice)
  {
  case 1:
    {
      imfrfList(Lat(L,1),reference_vertex,sequence);
    };break;
  case 2:
    {
      imfrfList(input_vertex,Lat(L,1),sequence);
    };break;
  case 3:
    {
      for (int i=1;i<=T1->getNbChild(input_vertex);i++)
      {
        int child_index = 1;
        int input_child = T1->child(input_vertex,i);
        if (!sameComplex(T1->getNode(input_child),T1->getNode(input_vertex)))
        {
          int reference_child = Lat(L,child_index++);
          if (reference_child != EMPTY_NODE)
          {
            TreeList(input_child,reference_child,sequence);
          }
        }
      }
    };break;
    default : break;
  }
}

DistanceType ExtendedMatching::match()
{
  DistanceType D=0.0;
  if ((!T1->isNull())&&(!T2->isNull()))
  {
    for (int input_vertex=T1->getNbVertex()-1;input_vertex>=0;input_vertex--)
    {
      _distances.openDistancesVector(input_vertex);

      for (int reference_vertex=T2->getNbVertex()-1;reference_vertex>=0;reference_vertex--)
      {
        distanceBetweenForest(input_vertex,reference_vertex);
        distanceBetweenTree(input_vertex,reference_vertex);
      }

      for (int i=1;i<=T1->getNbChild(input_vertex);i++)
      {
        _distances.closeDistancesVector(T1->child(input_vertex,i));
      }
    };

    D=getDBT(0,0);

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

DistanceType ExtendedMatching::getDBMF(int input_vertex,int reference_vertex) const
{
  return(_distances.getDBMF(input_vertex,reference_vertex));
}

DistanceType ExtendedMatching::getDBPF(int input_vertex,int reference_vertex) const
{
  return(_distances.getDBPF(input_vertex,reference_vertex));
}

DistanceType ExtendedMatching::getDBIMFRF(int input_vertex,int reference_vertex) const
{
  return(_distances.getDBIMFRF(input_vertex,reference_vertex));
}

DistanceType ExtendedMatching::getDBIFRMF(int input_vertex,int reference_vertex) const
{
  return(_distances.getDBIFRMF(input_vertex,reference_vertex));
}

DistanceType ExtendedMatching::getDBF(int input_vertex,int reference_vertex) const
{
  return(_distances.getDBF(input_vertex,reference_vertex));
}

DistanceType ExtendedMatching::getDBT(int input_vertex,int reference_vertex) const
{
  return(_distances.getDBT(input_vertex,reference_vertex));
}

int ExtendedMatching::M(int i_node,int r_node)
{
  return(_choices.getList(i_node,r_node)->back());
}

int ExtendedMatching::sameComplex(TreeNode* tree_node1,TreeNode* tree_node2)
{
  return(tree_node1->getComplex()==tree_node2->getComplex());
}







