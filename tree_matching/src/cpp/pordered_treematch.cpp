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


#include "pordered_treematch.h"

#include <math.h>
#include "tool/timer.h"
using namespace std;

TreeMatch_PO::TreeMatch_PO()
{
  _mtg                 = (MTG*) NULL;
  _roots               = (VIdList*) NULL;
  _localFun            = (NodeFunctionList*) NULL;
  _maxOrder            = 0;
  _maxNbVertex         = 0;
  _nbTree              = 0;
  _InsDelCostCoeff=.71;

}

TreeMatch_PO::TreeMatch_PO( MTG& mtg, Array* roots)
{
  _mtg=&mtg;
  assert(roots);
  assert((roots->surfaceType()==AMObjType::VTX)||(roots->surfaceType()==AMObjType::INTEGER));
  _roots=new VIdList();
  ArrayIter* pnext_root=roots->iterator();
  ArrayIter& next_root=*pnext_root;
  while(next_root())
    {
      if (next_root.key().tag() == AMObjType::VTX) _roots->insert(next_root.key().val.v);
      if (next_root.key().tag() == AMObjType::INTEGER) _roots->insert((VId) next_root.key().val.i);
    }

  _nbTree              = _roots->entries();

  _localFun            = (NodeFunctionList*) NULL;
  _InsDelCostCoeff     =.71;
  _distances.resize(_roots->entries());
  _time.resize(_roots->entries());
  _sequences.resize(_roots->entries());
  _sequences_tab.resize(_roots->entries());
  for (int tree=0;tree<_nbTree;tree++)
    {
      _distances[tree].resize(_nbTree-tree,0.0);
      _time[tree].resize(_nbTree-tree,0.0);
      _sequences[tree]=new SequenceVector(_nbTree-tree,(Sequence*) NULL);
      _sequences_tab[tree]=new SequenceTabVector(_nbTree-tree,(Sequence***) NULL);
    }

  topologicalMatching();

}



/* TreeMatching Topologique ordonné */

TreeMatch_PO::TreeMatch_PO(MTG& mtg,
			 Array* roots,
			 Array* local_functions,
			 AMString matching_type,
			 AMString mapping_type,
			 AMString mapping,
			 AMString scale_type,
			 const Vector_distance &ivect,
			 double coeff)
{
  _mtg=&mtg;
  // making of the vertex list
 _InsDelCostCoeff =coeff;

  _roots=new VIdList();

  assert(roots);
  assert((roots->surfaceType()==AMObjType::VTX)||(roots->surfaceType()==AMObjType::INTEGER));
  ArrayIter* pnext_root=roots->iterator();
  ArrayIter& next_root=*pnext_root;
  while(next_root())
    {
      if (next_root.key().tag() == AMObjType::VTX)
        {
          _roots->insert(next_root.key().val.v);
        }
      if (next_root.key().tag() == AMObjType::INTEGER)
        {
          _roots->insert((VId) next_root.key().val.i);
        }
    }

  _nbTree=_roots->entries();


  // making of the local functions list
  if (local_functions)
    {
      _vectorDist=*(new Vector_distance(ivect));
      _localFun=new NodeFunctionList();
      assert(local_functions->surfaceType()==AMObjType::FNODE);
      ArrayIter* p_next=local_functions->iterator();
      ArrayIter& next=*p_next;
      while(next())
        {
          if (next.key().tag() == AMObjType::FNODE)
            {
              AMModel* p_obj=next.key().val.p;
              if (p_obj)
                {
                  NodeFunction* Fun=new NodeFunction("Fun",(FNode*) next.key().val.p);
                  _localFun->push_back(*Fun);
                }
              else
                {
                  _localFun->push_back(*((NodeFunction*) NULL));
                }
            }
          else
            {
              _localFun->push_back(*((NodeFunction*)NULL));
            }
        }
    }
  else
    {
      _localFun = (NodeFunctionList*) NULL;
    }

  _distances.resize(_roots->entries());
  _time.resize(_roots->entries());
  
  _sequences.resize(_roots->entries());
  _sequences_tab.resize(_roots->entries());
  
  for (int tree=0;tree<_nbTree;tree++)
    {
      _distances[tree].resize(_nbTree-tree,0.0);
      _time[tree].resize(_nbTree-tree,0.0);
      _sequences[tree]=new SequenceVector(_nbTree-tree,(Sequence*) NULL);
      _sequences_tab[tree]=new SequenceTabVector(_nbTree-tree,(Sequence***) NULL);
    }

  

  if (!strcmp(matching_type,"Edition"))
    _matchingType = EDITION;
  if (!strcmp(matching_type,"Alignment"))
    _matchingType = ALIGNMENT;
  if (!strcmp(matching_type,"SmallestCommonSuperTree"))
    _matchingType = SMALLESTCOMMONSUPERTREE;
  if (!strcmp(matching_type,"LargestCommonSubTree"))
    _matchingType = LARGESTCOMMONSUBTREE;

  if (!strcmp(mapping_type,"Global"))
    _mappingType = GLOBAL;
  if (!strcmp(mapping_type,"Local"))
    _mappingType = LOCAL;

  if (!strcmp(mapping,"General"))
       _mapping = GENERAL;
  if (!strcmp(mapping,"Distance"))
    _mapping = DISTANCE;
  if (!strcmp(mapping,"Similarity"))
    _mapping = SIMILARITY;
  if (!strcmp(mapping,"EndSpaceFree"))
    _mapping = ENDSPACEFREE;
  
  

  if (!strcmp(scale_type,"SingleScale"))
    _scaleType = SINGLESCALE;
  if (!strcmp(mapping,"MultiScale"))
    _scaleType = MULTISCALE;

  if (_localFun)
    weightedMatching();
  else
    topologicalMatching();
}

//--------------------------------------------------------------------------------------
// DESTRUCTOR
//--------------------------------------------------------------------------------------


TreeMatch_PO::~TreeMatch_PO()
{
  if (_localFun)
    {
      _localFun->clear();
    }
  for (int it=0;it<_nbTree-1;it++)
    {
          if (_sequences[it])
            {
              for (int rt=0;rt<_nbTree-it;rt++)
                {
                  if ((*_sequences[it])[rt])
                    {
                      delete (Sequence*) (*_sequences[it])[rt];
                    }
                }
              delete (SequenceVector*) _sequences[it];
            }
    }

  for (int tree=0;tree<_nbTree;tree++)
    {
      if (_trees[tree]) delete (TreeGraph*) _trees[tree];
    }
}

//------------------------------------------------------------------------------------------------------
// TOPOLOGICAL MATCHING
//------------------------------------------------------------------------------------------------------
void TreeMatch_PO::topologicalMatching()
{
  int nb_tree=_roots->entries();
  _trees.resize(nb_tree);
  _maxNbVertex=0;
  _maxOrder=0;
  for (int i=0;i<nb_tree;i++)
    {
      VId root=_roots->at(i);
      _trees[i]=new TreeGraph(*_mtg,root,TOPO,ANY);
      _maxOrder = I_MAX(_maxOrder,_trees[i]->getOrder());
      _maxNbVertex=I_MAX(_trees[i]->getNbVertex(),_maxNbVertex);
    }
  for (int i_tree=0;i_tree<nb_tree;i_tree++)
    {
      for(int r_tree=i_tree+1;r_tree<nb_tree;r_tree++)
        {
          int tree_size2=_trees[r_tree]->getNbVertex();
          int tree_size1=_trees[i_tree]->getNbVertex();
          Sequence* matching_sequence=new Sequence();

          TOOLS(Timer) chrono;
          chrono.start();

          DistanceType matching_distance;

	    

	  matching_distance=MatchByTopology(*_trees[i_tree],
					    *_trees[r_tree],
					    matching_sequence);
	  double Time= chrono.stop();

	

	  putDistance(matching_distance,i_tree,r_tree);
	  putSequence(matching_sequence,i_tree,r_tree);
	  int seq_size=matching_sequence->getSize();
	  int sub_number=0;
	  int mat_number=0;
	  matching_sequence->reset();
	  do
	    {
	      if (matching_sequence->getCurrent()->getCost()==0.0)
		{
		  mat_number++;
		}
	      else
		{
		  sub_number++;
		}
	    } while(matching_sequence->next());
	  matching_sequence->putNbMat(mat_number);
	  matching_sequence->putNbSub(sub_number);
	  matching_sequence->putNbDel(tree_size1-seq_size);
	  matching_sequence->putNbIns(tree_size2-seq_size);
	  cout<<endl<<"MATCHING TIME "<<i_tree<<" with "<<r_tree<<" : "<<Time<<endl;
	  putTime(Time,i_tree,r_tree);
	  //Time= 0.0;
	}
    }
}


//------------------------------------------------------------------------------------------------------
// WEIGHTED MATCHING
//------------------------------------------------------------------------------------------------------
void TreeMatch_PO::weightedMatching()
{
  cerr<<"Weighted matching"<<endl;
  int nb_tree = _distances.size();
  int nb_fun  = _localFun->size();

  standardization();

  for (int i_tree=0;i_tree<nb_tree;i_tree++)
    {
      for(int r_tree=i_tree+1;r_tree<nb_tree;r_tree++)
        {
          Sequence* matching_sequence=new Sequence();

          TOOLS(Timer) chrono;
          chrono.start();

          DistanceType matching_distance=MatchByTopology(*_trees[i_tree],
                                                          *_trees[r_tree],
                                                          matching_sequence);
          double Time= chrono.stop();

          cout<<endl<<"MATCHING TIME "<<i_tree<<" with "<<r_tree<<" : "<<Time<<endl;
          putDistance(matching_distance,i_tree,r_tree);
          putSequence(matching_sequence,i_tree,r_tree);

          int seq_size=matching_sequence->getSize();
          int sub_number=0;
          int mat_number=0;
          matching_sequence->reset();
          do
            {
              if (matching_sequence->getCurrent()->getCost()<=0.01)
                {
                  mat_number++;
                }
              else
                {
                  sub_number++;
                }
            } while(matching_sequence->next());

          matching_sequence->putNbMat(mat_number);
          matching_sequence->putNbSub(sub_number);

          int tree_size2=_trees[r_tree]->getNbVertex();
          int tree_size1=_trees[i_tree]->getNbVertex();

          matching_sequence->putNbDel(tree_size1-seq_size);
          matching_sequence->putNbIns(tree_size2-seq_size);
        }
    }
}


//--------------------------------------------------------------------------------------
// DISTANCE FOR THE CLASSIC ALGORITHM
//--------------------------------------------------------------------------------------
DistanceType TreeMatch_PO::MatchByTopology(TreeGraph& Tree1,
					  TreeGraph& Tree2,
					  Sequence * S)
{

  NodeCost *MCF;
  Matching_PO *M;
  if (_localFun){
  if (_mappingType == LOCAL)
      MCF=new WeightedNodeCost(SCORE,_vectorDist,_dispersion,_maxValue,_minValue,_InsDelCostCoeff);
    else
    MCF = new WeightedNodeCost(WEIGTH,_vectorDist,_dispersion,_maxValue,_minValue,_InsDelCostCoeff);
  }
  else{
    if (_mappingType == LOCAL)
      MCF=new TopologicalLocalNodeCost(LOCAL_TOPO);
    else
      MCF=new NodeCost(TOPOLOGIC);
  }
  DistanceType D;
  if (_scaleType == SINGLESCALE){
    if (_matchingType == EDITION){
      if (_mappingType == GLOBAL){
	if (_mapping == GENERAL){
	  M = new Matching_PO(Tree1,Tree2,*MCF);
	}
	else {
	  if (_mapping == DISTANCE){
	    M = new Matching_PO(Tree1,Tree2,*MCF);
	  }
	  else {
	    if (_mapping == ENDSPACEFREE){
	      M = new Matching_PO_End_Space(Tree1,Tree2,*MCF);
	    }
	    else
	      cerr<<"Not yet Implemented"<<endl<<"General topological Matching"<<endl;
	  }
	}
      }
      
      else
	{
	  if (_mappingType == LOCAL){
	    cout<<"ici"<<endl;
	    M = new Matching_PO_Local(Tree1,Tree2,*MCF);
	  }
	  else
	    cerr<<"Not yet Implemented"<<endl<<"General topological Matching"<<endl;}
    }
    else 
      cerr<<"Not yet Implemented"<<endl<<"General topological Matching"<<endl;
  }
  else
    cerr<<"Not yet Implemented"<<endl<<"General topological Matching"<<endl;
  
  
  
  D = M->match();
  Sequence* s=new Sequence();
  
  M->TreeList(M->getI_v(),M->getR_v(),*s);

  int mat_number = 0;
  int sub_number = 0;
  DistanceType sub_cost = 0.0;
  //DistanceType del_cost = M->getDBT(0,EMPTY_TREE);
  DistanceType del_cost = 0;
  //DistanceType ins_cost = M->getDBT(EMPTY_TREE,0);
  DistanceType ins_cost = 0;
  if (s->getSize())
    {
      s->reset();
      do
        {
          TreeNode* tree_node1=Tree1.getNode(s->getCurrent()->getIV());
          TreeNode* tree_node2=Tree2.getNode(s->getCurrent()->getRV());
	  DistanceType scost = s->getCurrent()->getCost();
	  if (scost==0.0)
	    mat_number++;
	  else
	    sub_number++;
	  sub_cost += scost ;
	  del_cost -= MCF->getDeletionCost(tree_node1);
	  ins_cost -= MCF->getInsertionCost(tree_node2);
          S->append(tree_node1->getVertex(),tree_node2->getVertex(),MCF->getChangingCost(tree_node1,tree_node2));
        } while(s->next());
    }
  int ins_number = Tree2.getNbVertex() - (mat_number+sub_number);
  int del_number = Tree1.getNbVertex() - (mat_number+sub_number);

  S->putInsCost(ins_cost);
  S->putDelCost(del_cost);
  S->putSubCost(sub_cost);

  S->putNbIns(ins_number);
  S->putNbDel(del_number);
  S->putNbSub(sub_number);
 
  delete (Sequence*) s;
  delete (NodeCost*) MCF;
  delete (Matching_PO*) M;
  return(D);
}


//----------------------------------------
// STANDARTIZATION
//---------------------------------------
void TreeMatch_PO::standardization()
{
  //Methode de standardization
  // Si La norme est L1, on calcule l'abolute deviation (absd) pour chacune des fonctions
  //Elle correspond a la somme des differences des valeurs absolues pour toutes les valeurs
  //prise par les noueds
  //Sinon on calcule la standard deviation, qui remplace la valeur abolue par un carre
  int nb_tree = _distances.size();
  int nb_fun  = _localFun->size();

  DistanceType * attribute_range_max = new DistanceType [nb_fun];
  DistanceType * attribute_range_min = new DistanceType [nb_fun];

  _dispersion= ValueVector(nb_fun,MINDIST);
  _maxValue=  ValueVector(nb_fun,MINDIST);
  _minValue=  ValueVector(nb_fun,MINDIST);

  DistanceVector fun_norm(nb_fun,MINDIST);
  DistanceVector absd(nb_fun,MINDIST);


  int nbNode=0;
  _trees.resize(nb_tree);
  _maxNbVertex = 0;
  _maxOrder    = 0;
  int** frequency=new int*; // Nbre de valeur par variable
  TreeNode* tree_node=new TreeNode();
  int k,i;
  for (k=0;k<nb_fun;k++)
    {
      if (_vectorDist.get_variable_type(k)==SYMBOLIC)
        {
          frequency[k] = new int[_vectorDist.get_nb_value(k)];
          for ( i=0;i<_vectorDist.get_nb_value(k);i++)
            frequency[k][i]=0;
        }
      else
        {
          frequency[k]=NULL;
	  attribute_range_max[k] = MINDIST;
	  attribute_range_min[k] = MAXDIST;
        }
    }


  for ( i=0;i<nb_tree;i++)
    {
      VId root     = _roots->at(i);
      _trees[i]    = new TreeGraph(*_mtg,root,TOPO,ANY,TRUE,_localFun);
      _maxNbVertex = I_MAX(_maxNbVertex,_trees[i]->getNbVertex());
      _maxOrder    = I_MAX(_maxOrder,_trees[i]->getOrder());
      nbNode = nbNode + _trees[i]->getNbVertex();
      for (int j=0;j<_trees[i]->getNbVertex();j++)
        {
          tree_node=new TreeNode(*_trees[i]->getNode(j));
          for ( k=0;k<nb_fun;k++)
            {
              if (_vectorDist.get_variable_type(k)==SYMBOLIC)
                {
                  frequency[k][int(tree_node->getValue(k))]++;
                  //cout<<frequency[k][int(tree_node->getValue(k))]<<endl;
                }
              else
                {
                  DistanceType ccost =ABS(tree_node->getValue(k));
		  //cerr<<"Tree Node Value of "<<tree_node->getVertex()<<" = "<<ccost<<endl;
                  if (_vectorDist.get_distance_type()==ABSOLUTE_VALUE)
                    {
                      attribute_range_max[k]=MAX(attribute_range_max[k],ccost);
                      attribute_range_min[k]=MIN(attribute_range_min[k],ccost);
                    }
                  else
                    {
                      absd[k] = absd[k] + ccost*ccost;
                      fun_norm[k]=fun_norm[k]+ccost;
                      attribute_range_max[k]=MAX(attribute_range_max[k],ccost*ccost);
                      attribute_range_min[k]=MIN(attribute_range_min[k],ccost*ccost);
                    }
                }
            }
        }
    }
  for (k=0;k<nb_fun;k++)
    {
      if (_vectorDist.get_variable_type(k)==SYMBOLIC)
        {
          for (int i=0;i<_vectorDist.get_nb_value(k);i++)
            {
              absd[k]=absd[k]+double(frequency[k][i])*double(frequency[k][i])/nbNode/nbNode;
            }
          absd[k]=(1-absd[k])*nbNode/(nbNode-1);
        }
      else
        {
          if (_vectorDist.get_distance_type()==ABSOLUTE_VALUE)
            {
              absd[k] =  meanAbsoluteDifferenceComputation(k);
            }
          else
            {
              fun_norm[k]=fun_norm[k]*fun_norm[k];
              absd[k] = (2*(nbNode*absd[k]-fun_norm[k])/nbNode)/(nbNode-1);
            }
        }
      delete(frequency[k]);
    }
  for ( i=0;i<nb_tree;i++)
    {
      for (int j=0;j<_trees[i]->getNbVertex();j++)
        {
          //tree_node=_trees[i]->getNode(j);
          for ( k = 0;k<nb_fun;k++)
            {
              if (_vectorDist.get_variable_type(k)==SYMBOLIC)
                {
                  _dispersion[k]=absd[k];
                }
              else
                {
                  _dispersion[k]=absd[k];
                  _minValue[k]=attribute_range_min[k];
                  _maxValue[k]=attribute_range_max[k];
                }
            }
        }
      // _trees[i]->normalizeNodes(absd,indelcost[i]);
    }
  delete(frequency);
  delete(tree_node);
  delete attribute_range_min;
  delete attribute_range_max;
}



// Calul de la quantite de standardisation dans le cas valeur Absolue

DistanceType TreeMatch_PO::meanAbsoluteDifferenceComputation(int variable)
{
  int nb_tree = _distances.size();
  DistanceType mean_absolute_difference = 0;
  TreeNode* tree_node1=new TreeNode();
  TreeNode* tree_node2=new TreeNode();
  int nbNode=0;

  for (int i=0;i<nb_tree;i++)
    {
      nbNode+=_trees[i]->getNbVertex();
      for (int j1=0;j1<_trees[i]->getNbVertex();j1++)
        {
          tree_node1=new TreeNode(*_trees[i]->getNode(j1));
          for (int k1=j1+1;k1<_trees[i]->getNbVertex();k1++)
            {
              tree_node2=new TreeNode(*_trees[i]->getNode(k1));
              mean_absolute_difference += ABS(tree_node1->getValue(variable)-tree_node2->getValue(variable));
            }
        }
      for (int j=i+1;j<nb_tree;j++)
        {
          for (int k=0;k<_trees[i]->getNbVertex();k++)
            {
              tree_node1=new TreeNode(*_trees[i]->getNode(k));
              for (int m=0;m<_trees[j]->getNbVertex();m++)
                {
                  tree_node2=new TreeNode(*_trees[j]->getNode(m));
                  mean_absolute_difference += ABS(tree_node1->getValue(variable)-tree_node2->getValue(variable));
                }
            }
        }
    }

  delete(tree_node1);
  delete(tree_node2);
  mean_absolute_difference=2*mean_absolute_difference/(nbNode*(nbNode-1));
  return(mean_absolute_difference);

}
