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


#include "treematch.h"
#include "treenode.h"
#include <math.h>
#include "tool/timer.h"
using namespace std;
VPTOOLS_USING(Timer)

TreeMatch::TreeMatch()
{
  _mtg                 = (MTG*) NULL;
  _roots               = (VIdList*) NULL;
  _localFun            = (NodeFunctionList*) NULL;
  _maxOrder            = 0;
  _maxNbVertex         = 0;
  _nbTree              = 0;
  _InsDelCostCoeff=.71;

}

TreeMatch::TreeMatch( MTG& mtg, Array* roots)
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
  _intDistances.resize(_roots->entries());
  _sequences.resize(_roots->entries());
  _sequences_tab.resize(_roots->entries());
  for (int tree=0;tree<_nbTree;tree++)
    {
      _distances[tree].resize(_nbTree-tree,0.0);
      _time[tree].resize(_nbTree-tree,0.0);
      _intDistances[tree].resize(_nbTree-tree,0);
     _sequences[tree]=new SequenceVector(_nbTree-tree,(Sequence*) NULL);
      _sequences_tab[tree]=new SequenceTabVector(_nbTree-tree,(Sequence***) NULL);
    }

  topologicalMatching();

}


TreeMatch::TreeMatch(MTG& mtg,
                     Array* roots,
                     Array* local_functions,
                     AMString matching_type,
                     AMString ordered_type,
                     int self_similarity,
                     const Vector_distance &ivect,
                     double coeff)
{
  _scale = 3;
  _mtg=&mtg;
  // making of the vertex list
  _InsDelCostCoeff =coeff;

  _roots=new VIdList();
  _selfSimilarity = self_similarity;

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



  if (self_similarity != 1)
    {
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
    }

  if (!strcmp(matching_type,"by_weights"))
    _matchingType = BY_WEIGHTS;
  if (!strcmp(matching_type,"sequence"))
    _matchingType = SEQUENCE_MATCHING;
  if (!strcmp(matching_type,"selkow"))
    _matchingType = SELKOW_MATCHING;
  if (!strcmp(matching_type,"by_topology"))
    _matchingType = BY_TOPOLOGY;
  if (!strcmp(matching_type,"by_inclusion"))
    _matchingType = BY_INCLUSION;
  if (!strcmp(matching_type,"by_component"))
    _matchingType = BY_COMPONENTS;
  if (!strcmp(matching_type,"by_complex"))
    _matchingType = BY_COMPLEX;
  if (!strcmp(matching_type,"ordered"))
    _matchingType = ORDERED_MATCHING;
  if (!strcmp(matching_type,"test"))
    _matchingType = TEST;
  if (!strcmp(matching_type,"endSpaceFree"))
    _matchingType = END_SPACE_FREE;
  if (!strcmp(matching_type,"fractal"))
    _matchingType = FRACTAL;

  if (!strcmp(ordered_type,"jiangWangZhang"))
    _orderedType = JIANG_WANG_ZHANG;
  if (!strcmp(ordered_type,"TichitFerraro"))
    _orderedType = TICHIT_FERRARO;
 if (!strcmp(ordered_type,"TichitFerraro1"))
    _orderedType = TICHIT_FERRARO1;
 if (!strcmp(ordered_type,"FerraroOuangraoua"))
    _orderedType = FERRARO_OUANGRAOUA;
  if (!strcmp(ordered_type,"FerraroOuangraoua1"))
    _orderedType = FERRARO_OUANGRAOUA1;
  if (!strcmp(ordered_type,"FerraroOuangraoua2"))
    _orderedType = FERRARO_OUANGRAOUA2;
  if (!strcmp(ordered_type,"PartialOrder"))
    _orderedType = PARTIAL_ORDER;
  if (!strcmp(ordered_type,"ZhangShasha"))
    _orderedType = ZHANG_SHASHA;

  switch(_selfSimilarity)
    {
    case 0: {
      switch(_matchingType)
        {
        case BY_WEIGHTS       :{ weigthedMatching();    };break;
        case BY_TOPOLOGY      :{ topologicalMatching(); };break;
        case BY_INCLUSION     :{ topologicalMatching();   };break;
        case BY_COMPONENTS    :{ topologicalMatching();   };break;
        case BY_COMPLEX       :{ topologicalMatching();   };break;
        case ORDERED_MATCHING :{ topologicalMatching();   };break;
        case SEQUENCE_MATCHING :{ sequenceMatching();   };break;
        case SELKOW_MATCHING  :{ selkowMatching();   };break;
		case TEST             :{ topologicalMatching();   };break;
        default : assert(0);break;
        }
    }break;
    case 1: {
      switch(_matchingType)
        {
        case BY_WEIGHTS    :{ selfSimilarityWeigthedMatching();    };break;
        case BY_TOPOLOGY   :{ selfSimilarityTopologicalMatching(); };break;
        case BY_COMPLEX   :{ selfSimilarityComplexMatching(); };break;
        case END_SPACE_FREE   :{ endSpaceFreeTopologicalMatching(); };break;
        case FRACTAL   :{ fractalTopologicalMatching(); };break;
       default : assert(0);break;
        }
    }break;
    case 2: {
      switch(_matchingType)
        {
        case BY_WEIGHTS    :{ localWeigthedMatching();    };break;
        case BY_TOPOLOGY   :{ localTopologicalMatching(); };break;
        default : assert(0);break;
        }
    }break;
    default : assert(0);break;
    }

}


TreeMatch::TreeMatch(MTG& mtg,
                     Array* roots,
                     Array* local_functions,
                     AMString matching_type,
                     AMString ordered_type,
                     int self_similarity,
                     const Vector_distance &ivect,
                     double coeff,
                     NodeCost& Nd,
                     int scale)
{
  _scale = scale;
  _mtg=&mtg;
  // making of the vertex list
  _InsDelCostCoeff =coeff;
  _nd = &Nd;
  _roots=new VIdList();
  _selfSimilarity = self_similarity;
  
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



  if (self_similarity != 1)
    {
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
    }

  if (!strcmp(matching_type,"by_weights"))
    _matchingType = BY_WEIGHTS;
  if (!strcmp(matching_type,"sequence"))
    _matchingType = SEQUENCE_MATCHING;
  if (!strcmp(matching_type,"selkow"))
    _matchingType = SELKOW_MATCHING;
  if (!strcmp(matching_type,"by_topology"))
    _matchingType = BY_TOPOLOGY;
  if (!strcmp(matching_type,"by_inclusion"))
    _matchingType = BY_INCLUSION;
  if (!strcmp(matching_type,"by_component"))
    _matchingType = BY_COMPONENTS;
  if (!strcmp(matching_type,"by_complex"))
    _matchingType = BY_COMPLEX;
  if (!strcmp(matching_type,"ordered"))
    _matchingType = ORDERED_MATCHING;
  if (!strcmp(matching_type,"test"))
    _matchingType = TEST;
  if (!strcmp(matching_type,"endSpaceFree"))
    _matchingType = END_SPACE_FREE;
  if (!strcmp(matching_type,"fractal"))
    _matchingType = FRACTAL;

  if (!strcmp(ordered_type,"jiangWangZhang"))
    _orderedType = JIANG_WANG_ZHANG;
  if (!strcmp(ordered_type,"TichitFerraro"))
    _orderedType = TICHIT_FERRARO;
 if (!strcmp(ordered_type,"TichitFerraro1"))
    _orderedType = TICHIT_FERRARO1;
 if (!strcmp(ordered_type,"FerraroOuangraoua"))
    _orderedType = FERRARO_OUANGRAOUA;
  if (!strcmp(ordered_type,"FerraroOuangraoua1"))
    _orderedType = FERRARO_OUANGRAOUA1;
  if (!strcmp(ordered_type,"FerraroOuangraoua2"))
    _orderedType = FERRARO_OUANGRAOUA2;
  if (!strcmp(ordered_type,"PartialOrder"))
    _orderedType = PARTIAL_ORDER;
  if (!strcmp(ordered_type,"ZhangShasha"))
    _orderedType = ZHANG_SHASHA;

  switch(_selfSimilarity)
    {
    case 0: {
      switch(_matchingType)
        {
        case BY_WEIGHTS       :{ weigthedMatching();    };break;
        case BY_TOPOLOGY      :{ topologicalMatching(); };break;
        case BY_INCLUSION     :{ topologicalMatching();   };break;
        case BY_COMPONENTS    :{ topologicalMatching();   };break;
        case BY_COMPLEX       :{ topologicalMatching();   };break;
        case ORDERED_MATCHING :{ topologicalMatching();   };break;
        case SEQUENCE_MATCHING :{ sequenceMatching();   };break;
        case SELKOW_MATCHING  :{ selkowMatching();   };break;
		case TEST             :{ topologicalMatching();   };break;
        default : assert(0);break;
        }
    }break;
    case 1: {
      switch(_matchingType)
        {
        case BY_WEIGHTS    :{ selfSimilarityWeigthedMatching();    };break;
        case BY_TOPOLOGY   :{ selfSimilarityTopologicalMatching(); };break;
        case BY_COMPLEX   :{ selfSimilarityComplexMatching(); };break;
        case END_SPACE_FREE   :{ endSpaceFreeTopologicalMatching(); };break;
        case FRACTAL   :{ fractalTopologicalMatching(); };break;
       default : assert(0);break;
        }
    }break;
    case 2: {
      switch(_matchingType)
        {
        case BY_WEIGHTS    :{ localWeigthedMatching();    };break;
        case BY_TOPOLOGY   :{ localTopologicalMatching(); };break;
        default : assert(0);break;
        }
    }break;
    default : assert(0);break;
    }

}




//--------------------------------------------------------------------------------------
// DESTRUCTOR
//--------------------------------------------------------------------------------------


TreeMatch::~TreeMatch()
{
  if (_localFun)
    {
      _localFun->clear();
    }
  for (int it=0;it<_nbTree-1;it++)
    {
      if ( _matchingType != TEST)
        {
          if (_sequences[it])
            {
              /*for (int rt=0;rt<_nbTree-it;rt++)
                {
                  if ((*_sequences[it])[rt])
                    {
                      delete (Sequence*) (*_sequences[it])[rt];
                    }
		    }*/
              delete (SequenceVector*) _sequences[it];
            }
          // Il faut deleter _sequences_tab
        }
    }

  for (int tree=0;tree<_nbTree;tree++)
    {
      if (_trees[tree]) delete (TreeGraph*) _trees[tree];
    }
}


//--------------------------------------------------------------------------------------
// DISTANCE FOR THE INTERNAL ALGORITHM
//--------------------------------------------------------------------------------------
DistanceType TreeMatch::MatchByAttribute(TreeGraph& Tree1,
                                         TreeGraph& Tree2,
                                         Sequence * S)

{


  // making of the appropriate node function
  NodeCost* MCF;
  

  MCF=new WeightedNodeCost(WEIGTH,_vectorDist,_dispersion,_maxValue,_minValue,_InsDelCostCoeff);

  // making of the matching
  Matching* M=new Matching(Tree1,Tree2,*MCF);

  // effective matching : the result is stored in D
  DistanceType D=M->match();
  // making of the matching list
  Sequence s=Sequence();
  M->TreeList(0,0,s);
  //  cout<<"Construction de la séquence"<<endl;
  S->putInsCost(M->getInsertCost());
  S->putDelCost(M->getDeleteCost());
  S->putSubCost(M->getSubstitutionCost());
  if (s.getSize())
    {
      s.reset();
      do
        {
              TreeNode* tree_node1=Tree1.getNode(s.getCurrent()->getIV());
              TreeNode* tree_node2=Tree2.getNode(s.getCurrent()->getRV());
              S->append(tree_node1->getVertex(),tree_node2->getVertex(),((WeightedNodeCost*)MCF)->getChangingCost(tree_node1,tree_node2));
        } while(s.next());
    }
  //  delete (Sequence*) s;
  delete (NodeCost*) MCF;
  delete (Matching*) M;
  return(D);
}
//--------------------------------------------------------------------------------------
// LOCAL DISTANCE 
//--------------------------------------------------------------------------------------
DistanceType TreeMatch::LocalMatchByAttribute(TreeGraph& Tree1,
                                         TreeGraph& Tree2,
                                         Sequence * S)

{


  // making of the appropriate node function
  NodeCost* MCF;

  MCF=new WeightedNodeCost(WEIGTH,_vectorDist,_dispersion,_maxValue,_minValue,_InsDelCostCoeff);

  // making of the matching
  LocalMatching* M=new LocalMatching(Tree1,Tree2,*MCF);

  // effective matching : the result is stored in D
  DistanceType D=M->match();
  // making of the matching list
  Sequence* s=new Sequence();
  M->TreeList(0,0,*s);
  if (s->getSize())
    {
      s->reset();
      do
        {
              TreeNode* tree_node1=Tree1.getNode(s->getCurrent()->getIV());
              TreeNode* tree_node2=Tree2.getNode(s->getCurrent()->getRV());
              S->append(tree_node1->getVertex(),tree_node2->getVertex(),((WeightedNodeCost*)MCF)->getChangingCost(tree_node1,tree_node2));
        } while(s->next());
    }
  delete (Sequence*) s;
  delete (NodeCost*) MCF;
  delete (LocalMatching*) M;
  return(D);
}
//--------------------------------------------------------------------------------------
// DISTANCE FOR THE INTERNAL ALGORITHM
//--------------------------------------------------------------------------------------
DistanceType TreeMatch::MatchByAttribute(TreeGraph& Tree1,
                                         TreeGraph& Tree2)
{

  cout<<"Match by Attribute"<<endl;
  // making of the appropriate node function
  NodeCost* MCF;

  int nb_fun  = _localFun->size();
 _dispersion= ValueVector(nb_fun,1.);
  _maxValue=  ValueVector(nb_fun,0.);
  _minValue=  ValueVector(nb_fun,0.);

  MCF=new WeightedNodeCost(WEIGTH,_vectorDist,_dispersion,_maxValue,_minValue,_InsDelCostCoeff);

  // making of the matching
  Matching* M=new Matching(Tree1,Tree2,*MCF);
  // effective matching : the result is stored in D
  DistanceType D=M->match();

  for (int i = 0; i<Tree1.getNbVertex();i++)
    {
       //TreeNode* tree_node = Tree1.getNode(i);
      if ((Tree1.getNode(i))->getOrder() == 0)
	_trunk.push_back(i);
      //delete(tree_node);
     for (int j = i+1; j<Tree1.getNbVertex();j++)
        {
          DistanceType matching_distance = M->getDistanceMatrix(i,j);
          putDistance(matching_distance,i,j);
          /*Sequence* matching_sequence=new Sequence();
          Sequence* s = new Sequence();
          M->TreeList(i,j,*s);

          if (s->getSize())
            {
              s->reset();
              do
                {
                  TreeNode* tree_node1=Tree1.getNode(s->getCurrent()->getIV());
                  TreeNode* tree_node2=Tree2.getNode(s->getCurrent()->getRV());
                  matching_sequence->append(tree_node1->getVertex(),tree_node2->getVertex(),MCF->getChangingCost(tree_node1,tree_node2));
                } while(s->next());
            }
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

          int tree_size2=Tree1.getNbDesc(j) +1;
          int tree_size1=Tree1.getNbDesc(i) +1;

          matching_sequence->putNbDel(tree_size1-seq_size);
          matching_sequence->putNbIns(tree_size2-seq_size);
          putSequence(matching_sequence,i,j);

          // delete (Sequence*) matching_sequence;
          delete (Sequence*) s;*/
        }
    }
  delete (NodeCost*) MCF;
  //delete (Matching*) M;
  return(D);
}

//--------------------------------------------------------------------------------------
// DISTANCE Test
//--------------------------------------------------------------------------------------


DistanceType TreeMatch::MatchTest(TreeGraph& Tree1,
                                  TreeGraph& Tree2,
                                  Sequence * S)

{


  // making of the appropriate node function
  NodeCost* MCF;
  // making of the matching
  MCF=new NodeCost(TOPOLOGIC);

  MultiscaleMatching* M=new MultiscaleMatching(Tree1,Tree2,*MCF);


  // effective matching : the result is stored in D
  DistanceType D=M->match();
  // making of the matching list
  Sequence* s=new Sequence();
  M->TreeList(0,0,*s);
  if (s->getSize())
    {
      s->reset();
      do
        {
          TreeNode* tree_node1=Tree1.getNode(s->getCurrent()->getIV());
          TreeNode* tree_node2=Tree2.getNode(s->getCurrent()->getRV());
          S->append(tree_node1->getVertex(),tree_node2->getVertex(),MCF->getChangingCost(tree_node1,tree_node2));
        } while(s->next());
    }
  delete (Sequence*) s;
  delete (NodeCost*) MCF;
  delete (MultiscaleMatching*) M;
  return(D);
}


//--------------------------------------------------------------------------------------
// Inclusion d'un arbre dans un autre
//--------------------------------------------------------------------------------------
DistanceType TreeMatch::MatchByInclusion( TreeGraph& Tree1,
                                          TreeGraph& Tree2,
                                          Sequence * S)
{

  // NodeCostFunction used
  NodeCost* MCF=new NodeCost(TOPOLOGIC);
  // making of the matching
  MatchingMinimizeComponents* M=new MatchingMinimizeComponents(Tree1,Tree2,*MCF);
  // effective matching : the result is stored in D
  DistanceType D=M->match();
  // making of the matching list
  Sequence* s=new Sequence();
  M->TreeList(0,0,*s);
  //  S->putInsCost(M->getInsertCost());
  //  S->putDelCost(M->getDeleteCost());
  //  S->putSubCost(M->getSubstitutionCost());

  if (s->getSize())
    {
      s->reset();
      do
        {
          TreeNode* tree_node1=Tree1.getNode(s->getCurrent()->getIV());
          TreeNode* tree_node2=Tree2.getNode(s->getCurrent()->getRV());
          S->append(tree_node1->getVertex(),tree_node2->getVertex(),MCF->getChangingCost(tree_node1,tree_node2));
        } while(s->next());
    }
  delete (NodeCost*) MCF;
  delete (Sequence*) s;
  delete (MatchingMinimizeComponents*) M;
  return(D);
}
//--------------------------------------------------------------------------------------
// Ordered matching Jiang Zhang Algorithm
//--------------------------------------------------------------------------------------
DistanceType TreeMatch::JiangWangZhangMatching(TreeGraph& Tree1,
                                        TreeGraph& Tree2,
                                        Sequence * S)
{

  // NodeCostFunction used
  NodeCost* MCF=new NodeCost(TOPOLOGIC);
  JiangWangZhangMetric2* M=new JiangWangZhangMetric2(Tree1,Tree2);
  //effective matching : the result is stored in D
  DistanceType D=M->run();
  Sequence* s=new Sequence();
  s = M->getSequence();
//   S->putInsCost(M->getInsertCost());
//   S->putDelCost(M->getDeleteCost());
//   S->putSubCost(M->getSubstitutionCost());

  if (s->getSize())
    {
      s->reset();
      do
        {
          TreeNode* tree_node1=Tree1.getNode(s->getCurrent()->getIV());
          TreeNode* tree_node2=Tree2.getNode(s->getCurrent()->getRV());
          S->append(tree_node1->getVertex(),tree_node2->getVertex(),MCF->getChangingCost(tree_node1,tree_node2));
        } while(s->next());
    }
  delete (NodeCost*) MCF;
  delete (Sequence*) s;
  delete (JiangWangZhangMetric2*) M;
  return(-1.0*D);
}

//--------------------------------------------------------------------------------------
// Ordered matching Jiang Zhang Algorithm
//--------------------------------------------------------------------------------------
DistanceType TreeMatch::TichitFerraroMatching(TreeGraph& Tree1,
                                        TreeGraph& Tree2,
                                        Sequence * S,int scale)
{

  // NodeCostFunction used
  NodeCost* MCF=new NodeCost(TOPOLOGIC);
  TichitFerraro* M=new TichitFerraro(*_mtg,Tree1,Tree2,scale,1);
  //effective matching : the result is stored in D
  DistanceType D=M->run();
  Sequence* s=new Sequence();
  s = M->getSequence();
//   S->putInsCost(M->getInsertCost());
//   S->putDelCost(M->getDeleteCost());
//   S->putSubCost(M->getSubstitutionCost());

  if (s->getSize())
    {
      s->reset();
      do
        {
          TreeNode* tree_node1=Tree1.getNode(s->getCurrent()->getIV());
          TreeNode* tree_node2=Tree2.getNode(s->getCurrent()->getRV());
          S->append(tree_node1->getVertex(),tree_node2->getVertex(),MCF->getChangingCost(tree_node1,tree_node2));
        } while(s->next());
    }
  delete (NodeCost*) MCF;
  delete (Sequence*) s;
  delete (TichitFerraro*) M;
  return(-1.0*D);
}

DistanceType TreeMatch::TichitFerraro1Matching(TreeGraph& Tree1,
                                        TreeGraph& Tree2,
                                        Sequence * S,int scale)
{

  // NodeCostFunction used
  NodeCost* MCF=new NodeCost(TOPOLOGIC);
  TichitFerraro1* M=new TichitFerraro1(*_mtg,Tree1,Tree2,scale,1);
  //effective matching : the result is stored in D
  DistanceType D=M->run();
  Sequence* s=new Sequence();
  s = M->getSequence();
//   S->putInsCost(M->getInsertCost());
//   S->putDelCost(M->getDeleteCost());
//   S->putSubCost(M->getSubstitutionCost());

  if (s->getSize())
    {
      s->reset();
      do
        {
          TreeNode* tree_node1=Tree1.getNode(s->getCurrent()->getIV());
          TreeNode* tree_node2=Tree2.getNode(s->getCurrent()->getRV());
          S->append(tree_node1->getVertex(),tree_node2->getVertex(),MCF->getChangingCost(tree_node1,tree_node2));
        } while(s->next());
    }
  delete (NodeCost*) MCF;
  delete (Sequence*) s;
  delete (TichitFerraro1*) M;
  return(-1.0*D);
}

//--------------------------------------------------------------------------------------
// Ordered matching Jiang Zhang Multiscale Algorithm
//--------------------------------------------------------------------------------------

DistanceType TreeMatch::FerraroOuangraouaMatching(TreeGraph& Tree1,
                                        TreeGraph& Tree2,
                                        Sequence * S,int scale)
{
  // NodeCostFunction used
  NodeCost* MCF=new NodeCost(TOPOLOGIC);
  //FerraroOuangraoua* M=new FerraroOuangraoua(Tree1,Tree2,*_nd, scale);
  FerraroOuangraoua* M=new FerraroOuangraoua(Tree1,Tree2,*MCF, scale);
  //effective matching : the result is stored in D
  DistanceType D=M->run();
  Sequence* s = M->getSequence();
  S->putInsCost(M->getInsertionCost());
  S->putDelCost(M->getDeletionCost());
  S->putNbIns(M->getSequence()->getNbIns());
  S->putNbDel(M->getSequence()->getNbDel());
  S->putSubCost(M->getMatchingCost());
  //cout<<"OK1 "<<endl;
  if (s->getSize())
    {
      s->reset();
      do
	{
	  S->append(s->getCurrent()->getIV(),s->getCurrent()->getRV(),s->getCurrent()->getCost());
	} while(s->next());
      
    } 
  //cout<<"OK2 "<<endl;
  delete (NodeCost*) MCF;
  
  delete (FerraroOuangraoua*) M;
  return(D);
}

DistanceType TreeMatch::FerraroOuangraoua1Matching(TreeGraph& Tree1,
                                        TreeGraph& Tree2,
                                        Sequence * S,int scale)
{
  
  // NodeCostFunction used
  NodeCost* MCF=new NodeCost(TOPOLOGIC);
  FerraroOuangraoua1* M=new FerraroOuangraoua1(Tree1,Tree2,*MCF, scale);
  //effective matching : the result is stored in D
  DistanceType D=M->run();

  Sequence* s = M->getSequence();
  S->putInsCost(M->getInsertionCost());
  S->putDelCost(M->getDeletionCost());
  S->putSubCost(M->getMatchingCost());
  //cout<<"OK1 "<<endl;
  if (s->getSize())
    {
      s->reset();
      do
	{
	  S->append(s->getCurrent()->getIV(),s->getCurrent()->getRV(),s->getCurrent()->getCost());
	} while(s->next());
      
    } 
  //cout<<"OK2 "<<endl;
  delete (NodeCost*) MCF;
  
  delete (FerraroOuangraoua1*) M;
  return(D);
}

DistanceType TreeMatch::FerraroOuangraoua2Matching(TreeGraph& Tree1,
                                        TreeGraph& Tree2,
                                        Sequence * S,int scale)
{
  
  // NodeCostFunction used
  NodeCost* MCF=new NodeCost(TOPOLOGIC);
  FerraroOuangraoua2* M=new FerraroOuangraoua2(Tree1,Tree2,*MCF, scale);
  //effective matching : the result is stored in D
  DistanceType D=M->run();

  Sequence* s = M->getSequence();
  S->putInsCost(M->getInsertionCost());
  S->putDelCost(M->getDeletionCost());
  S->putSubCost(M->getMatchingCost());
  //cout<<"OK1 "<<endl;
  if (s->getSize())
    {
      s->reset();
      do
	{
	  S->append(s->getCurrent()->getIV(),s->getCurrent()->getRV(),s->getCurrent()->getCost());
	} while(s->next());
      
    } 
  //cout<<"OK2 "<<endl;
  delete (NodeCost*) MCF;
  
  delete (FerraroOuangraoua2*) M;
  return(D);
}

DistanceType TreeMatch::MatchPartialOrder(TreeGraph& Tree1,
					  TreeGraph& Tree2,
					  Sequence * S)
{
  
  
  // NodeCostFunction used
  NodeCost* MCF=new NodeCost(TOPOLOGIC);
  //  int optimization = 1;
  MatchingGposde* M = new MatchingGposde(Tree1,Tree2,*MCF);
  //effective matching : the result is stored in D
  
  DistanceType D=M->match();

 
  Sequence* s=new Sequence();
  M->TreeList(0,0,*s);
  S->putInsCost(M->getInsertionCost());
  S->putDelCost(M->getDeletionCost());
  S->putSubCost(M->getMatchingCost());
  cout<<"Sequence :"<<endl;
  if (s->getSize())
    {
      s->reset();
      do
        {
          TreeNode* tree_node1=Tree1.getNode(s->getCurrent()->getIV());
          TreeNode* tree_node2=Tree2.getNode(s->getCurrent()->getRV());
          S->append(tree_node1->getVertex(),tree_node2->getVertex(),MCF->getChangingCost(tree_node1,tree_node2));
	  cout<<s->getCurrent()->getIV()<<";"<<s->getCurrent()->getRV()<<endl;
        } while(s->next());
    }
  delete (NodeCost*) MCF;
  delete (MatchingGposde*) M;
  //  cout<<"End Match by Topology "<<endl;
  return(D);
}


//--------------------------------------------------------------------------------------
// DISTANCE FOR THE CLASSIC ALGORITHM
//--------------------------------------------------------------------------------------
DistanceType TreeMatch::MatchByTopology(TreeGraph& Tree1,
                                        TreeGraph& Tree2,
                                        Sequence * S)
{

  //  cout<<" Match by Topology "<<endl;
  // NodeCostFunction used
  NodeCost* MCF=new NodeCost(TOPOLOGIC);
  //  int optimization = 1;
  Matching* M=new Matching(Tree1,Tree2,*MCF);
  //effective matching : the result is stored in D
  DistanceType D=M->match();
  Sequence* s=new Sequence();
  M->TreeList(0,0,*s);
  S->putInsCost(M->getInsertCost());
  S->putDelCost(M->getDeleteCost());
  S->putSubCost(M->getSubstitutionCost());

  if (s->getSize())
    {
      s->reset();
      do
        {
          TreeNode* tree_node1=Tree1.getNode(s->getCurrent()->getIV());
          TreeNode* tree_node2=Tree2.getNode(s->getCurrent()->getRV());
          S->append(tree_node1->getVertex(),tree_node2->getVertex(),MCF->getChangingCost(tree_node1,tree_node2));
        } while(s->next());
    }
  delete (NodeCost*) MCF;
  delete (Sequence*) s;
  delete (Matching*) M;
  //  cout<<"End Match by Topology "<<endl;
  return(D);
}
//--------------------------------------------------------------------------------------
// DISTANCE FOR THE LOCAL ALGORITHM
//--------------------------------------------------------------------------------------
DistanceType TreeMatch::LocalMatchByTopology(TreeGraph& Tree1,
                                        TreeGraph& Tree2,
                                        Sequence * S)
{
  // NodeCostFunction used
  NodeCost* MCF=new TopologicalLocalNodeCost(LOCAL_TOPO);
  LocalMatching* M=new LocalMatching(Tree1,Tree2,*MCF);
  //effective matching : the result is stored in D
  DistanceType D=M->match();
  Sequence* s=new Sequence();
  M->TreeList(M->getMaxInput(),M->getMaxRef(),*s);

  if (s->getSize())
    {
      s->reset();
      do
        {
          TreeNode* tree_node1=Tree1.getNode(s->getCurrent()->getIV());
          TreeNode* tree_node2=Tree2.getNode(s->getCurrent()->getRV());
          S->append(tree_node1->getVertex(),tree_node2->getVertex(),MCF->getChangingCost(tree_node1,tree_node2));
        } while(s->next());
    }
  delete (NodeCost*) MCF;
  delete (Sequence*) s;
  delete (LocalMatching*) M;
  return(D);
}
//--------------------------------------------------------------------------------------
// Inclusion d'un arbre dans un autre
//--------------------------------------------------------------------------------------
DistanceType TreeMatch::MatchByTopology( TreeGraph& Tree1,
                                         TreeGraph& Tree2)
{

  cout<<"Match by Topology "<<endl;
  // NodeCostFunction used
  NodeCost* MCF=new NodeCost(TOPOLOGIC);
  // making of the matching
  int optimization = 1;
  //Matching* M=new Matching(Tree1,Tree2,*MCF ,optimization);
  SelfSimilarMatching* M=new SelfSimilarMatching(Tree1,Tree2,*MCF );
  // effective matching : the result is stored in D
  DistanceType D=M->match();
  int size1 = Tree1.getNbVertex();
  int size2 = Tree2.getNbVertex();
  // making of the matching list
//  cerr << "\x0d" <<"Already Saved : 0% ...   " << flush;
  int nb_tree = _roots->entries();
  if (nb_tree == 1){
    for (int i = 0; i<size1;i++)
      {
	//TreeNode* tree_node = Tree1.getNode(i);
	if ((Tree1.getNode(i))->getOrder() == 0)
	  _trunk.push_back(i);
	//delete(tree_node);
	for (int j = i+1 ; j<size2;j++)
	  {
	    DistanceType matching_distance = M->getDistanceMatrix(i,j);
	    putDistance(matching_distance,i,j);
	  }
	//if (int(100.*i/size1)%5 == 0)
	  //cerr << "\x0d" <<"Already Saved : "<< int(100.*i/size1) <<"% " <<" ...  " << flush;
      } 
  }
  else{
    for (int i = 0; i<size1;i++)
      {
	//TreeNode* tree_node = Tree1.getNode(i);
	if ((Tree1.getNode(i))->getOrder() == 1)
	  _trunk.push_back(i);
	//delete(tree_node);
	for (int j = 0 ; j<size2;j++)
	  {
	    DistanceType matching_distance = M->getDBT(i,j);
	    putDistance(matching_distance,i,j+size1);
	  }
	if (int(100.*i/size1)%5 == 0)
	  cerr << "\x0d" <<"Already Saved : "<< int(100.*i/size1) <<"% " <<" ...  " << flush;
      } 
  }
  cerr<<"\x0d"<<flush;
  delete (NodeCost*) MCF;
  delete (SelfSimilarMatching*) M;
  return(D);
}

//--------------------------------------------------------------------------------------
// Matching Using Complex [fer2000]
//--------------------------------------------------------------------------------------

DistanceType TreeMatch::MatchByComplex(TreeGraph& Tree1,
                                       TreeGraph& Tree2,
                                       Sequence * S)
{

  cout<<" Match by Complex "<<endl;
  // NodeCostFunction used
  NodeCost* MCF=new NodeCost(TOPOLOGIC);
  // making of the matching
  MatchingWithComplex* M=new MatchingWithComplex(Tree1,Tree2,*MCF);
  // effective matching : the result is stored in D
  DistanceType D=M->match();
  // making of the matching list
  Sequence* s=new Sequence();
  M->TreeList(0,0,*s);

  if (s->getSize())
    {
      s->reset();
      do
        {
          TreeNode* tree_node1=Tree1.getNode(s->getCurrent()->getIV());
          TreeNode* tree_node2=Tree2.getNode(s->getCurrent()->getRV());
          S->append(tree_node1->getVertex(),tree_node2->getVertex(),MCF->getChangingCost(tree_node1,tree_node2));
        } while(s->next());
    }



  delete (NodeCost*) MCF;
  delete (Sequence*) s;
  delete (MatchingWithComplex*) M;
  return(D);
}
//--------------------------------------------------------------------------------------
// Inclusion d'un arbre dans un autre
//--------------------------------------------------------------------------------------
DistanceType TreeMatch::MatchByComplex( TreeGraph& Tree1,
                                         TreeGraph& Tree2)
{

  cout<<"Match by Complex "<<endl;
  // NodeCostFunction used
  NodeCost* MCF=new NodeCost(TOPOLOGIC);
  // making of the matching
  int optimization = 1;
  //Matching* M=new Matching(Tree1,Tree2,*MCF ,optimization);
  SelfSimilarComplexMatching* M=new SelfSimilarComplexMatching(Tree1,Tree2,*MCF );
  // effective matching : the result is stored in D
  DistanceType D=M->match();
  int size1 = Tree1.getNbVertex();
  int size2 = Tree2.getNbVertex();
  // making of the matching list
//  cerr << "\x0d" <<"Already Saved : 0% ...   " << flush;
  int nb_tree = _roots->entries();
  if (nb_tree == 1){
    for (int i = 0; i<size1;i++)
      {
	//TreeNode* tree_node = Tree1.getNode(i);
		if ((Tree1.getNode(i))->getOrder() == 0)
			_trunk.push_back(i);
	//delete(tree_node);
    	for (int j = 0 ; j<size2;j++)
		{
		    DistanceType matching_distance = M->getDistanceMatrix(i,j);
		    putDistance(matching_distance,i,j);
		}
	//if (int(100.*i/size1)%5 == 0)
	  //cerr << "\x0d" <<"Already Saved : "<< int(100.*i/size1) <<"% " <<" ...  " << flush;
      } 
  }
  else{
    for (int i = 0; i<size1;i++)
      {
	//TreeNode* tree_node = Tree1.getNode(i);
	if ((Tree1.getNode(i))->getOrder() == 1)
	  _trunk.push_back(i);
	//delete(tree_node);
	for (int j = 0 ; j<size2;j++)
	  {
	    DistanceType matching_distance = M->getDBT(i,j);
	    putDistance(matching_distance,i,j+size1);
	  }
	if (int(100.*i/size1)%5 == 0)
	  cerr << "\x0d" <<"Already Saved : "<< int(100.*i/size1) <<"% " <<" ...  " << flush;
      } 
  }
  cerr<<"\x0d"<<flush;
  delete (NodeCost*) MCF;
  delete (SelfSimilarComplexMatching*) M;
  return(D);
}

// ------------------------------------------------------
// Matching Minimizing the number of connected components
//-------------------------------------------------------


DistanceType TreeMatch::MatchByComponents(TreeGraph& Tree1,
                                          TreeGraph& Tree2,
                                          Sequence * S)
{

  cout<<" Match by Connected Components "<<endl;
  // NodeCostFunction used
  NodeCost* MCF=new NodeCost(TOPOLOGIC);
  // making of the matching
  MatchingByComponents* M=new MatchingByComponents(Tree1,Tree2,*MCF);
  // effective matching : the result is stored in D
  DistanceType D=M->match();
  // cout<<"End Matching "<<endl;


  // making of the matching list
  Sequence* s=new Sequence();
  M->TreeList(0,0,*s);

  if (s->getSize())
    {
      s->reset();
      do
        {
          TreeNode* tree_node1=Tree1.getNode(s->getCurrent()->getIV());
          TreeNode* tree_node2=Tree2.getNode(s->getCurrent()->getRV());
          S->append(tree_node1->getVertex(),tree_node2->getVertex(),MCF->getChangingCost(tree_node1,tree_node2));
        } while(s->next());
    }
  S->putInsCost(M->getInsertCost());
  S->putDelCost(M->getDeleteCost());
  S->putSubCost(M->getSubstitutionCost());

  delete (NodeCost*) MCF;
  delete (Sequence*) s;
  delete (Matching*) M;
  return(D);
}


//------------------------------------------------------------------------------------------------------
// WEIGHTED MATCHING
//------------------------------------------------------------------------------------------------------
void TreeMatch::weigthedMatching()
{
  int nb_tree = _distances.size();
  int nb_fun  = _localFun->size();

  standardization();

  for (int i_tree=0;i_tree<nb_tree;i_tree++)
    {
      for(int r_tree=i_tree+1;r_tree<nb_tree;r_tree++)
        {
          Sequence* matching_sequence=new Sequence();

          Timer chrono;
          chrono.start();

          DistanceType matching_distance=MatchByAttribute(*_trees[i_tree],
                                                          *_trees[r_tree],
                                                          matching_sequence);
          double Time= chrono.stop();

          cout<<"MATCHING TIME "<<i_tree<<" with "<<r_tree<<" : "<<Time<<endl;
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

//------------------------------------------------------------------------------------------------------
// SELF SIMILAR WEIGHTED MATCHING
//------------------------------------------------------------------------------------------------------
void TreeMatch::selfSimilarityWeigthedMatching()
{
  int nb_fun  = _localFun->size();
  int nb_tree=_roots->entries();
  _trees.resize(nb_tree);
  _maxNbVertex=0;
  _maxOrder=0;
  for (int i=0;i<nb_tree;i++)
    {
      VId root=_roots->at(i);
      _trees[i]=new TreeGraph(*_mtg,root,TOPO,ANY,TRUE,_localFun);
      _maxOrder = I_MAX(_maxOrder,_trees[i]->getOrder());
      _maxNbVertex=I_MAX(_trees[i]->getNbVertex(),_maxNbVertex);
    }

//   cout<<"standardization"<<endl;
//   standardization();
//   cout<<"end standardization"<<endl;

  int nbTree = _trees[0]->getNbVertex();
  _time.resize(2);
  _distances.resize(nbTree);
  _sequences.resize(nbTree);
  for (int tree=0;tree<nbTree;tree++)
    {
      //cout<<tree<<endl;
       _sequences[tree]=new SequenceVector(nbTree-tree,(Sequence*) NULL);
       //cout<<tree<<endl;
      _distances[tree].resize(nbTree-tree,0.0);
    }

  Timer chrono;
  chrono.start();

  DistanceType matching_distance=MatchByAttribute(*_trees[0],
                                                  *_trees[0]);
  double Time = chrono.stop();
  cout<<"MATCHING TIME : "<<Time<<endl;
}

//------------------------------------------------------------------------------------------------------
// LOCAL WEIGHTED MATCHING
//------------------------------------------------------------------------------------------------------
void TreeMatch::localWeigthedMatching()
{
  int nb_tree = _distances.size();
  int nb_fun  = _localFun->size();

  standardization();

  for (int i_tree=0;i_tree<nb_tree;i_tree++)
    {
      for(int r_tree=i_tree+1;r_tree<nb_tree;r_tree++)
        {
          Sequence* matching_sequence=new Sequence();

          Timer chrono;
          chrono.start();

          DistanceType matching_distance=LocalMatchByAttribute(*_trees[i_tree],
                                                          *_trees[r_tree],
                                                          matching_sequence);
          double Time= chrono.stop();

          cout<<"MATCHING TIME "<<i_tree<<" with "<<r_tree<<" : "<<Time<<endl;
          putDistance(matching_distance,i_tree,r_tree);
          putSequence(matching_sequence,i_tree,r_tree);
		  putTime(Time,i_tree,r_tree);

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



//----------------------------------------
// STANDARTIZATION
//---------------------------------------
void TreeMatch::standardization()
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

DistanceType TreeMatch::meanAbsoluteDifferenceComputation(int variable)
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

//------------------------------------------------------------------------------------------------------
// TOPOLOGICAL MATCHING
//------------------------------------------------------------------------------------------------------
void TreeMatch::topologicalMatching()
{
  //cout<<"topological matching"<<endl;
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

          Timer chrono;
          chrono.start();

          DistanceType matching_distance;
          if (_matchingType==BY_TOPOLOGY)
            {
              matching_distance=MatchByTopology(*_trees[i_tree],
                                                *_trees[r_tree],
                                                matching_sequence);
            }
          if (_matchingType==BY_INCLUSION)
            {
              matching_distance=MatchByInclusion(*_trees[i_tree],
                                                 *_trees[r_tree],
                                                 matching_sequence);
            }
          if (_matchingType==BY_COMPONENTS)
            {
              matching_distance=MatchByComponents(*_trees[i_tree],
                                                  *_trees[r_tree],
                                                  matching_sequence);
            }
          if (_matchingType==BY_COMPLEX)
            {
              matching_distance=MatchByComplex(*_trees[i_tree],
                                               *_trees[r_tree],
                                               matching_sequence);
            }
	  if (_matchingType==ORDERED_MATCHING)
	    {
	      if (_orderedType == JIANG_WANG_ZHANG)
		matching_distance = JiangWangZhangMatching(*_trees[i_tree],
							    *_trees[r_tree],
							    matching_sequence);
	      else if (_orderedType == TICHIT_FERRARO)
		matching_distance = TichitFerraroMatching(*_trees[i_tree],
							    *_trees[r_tree],
							    matching_sequence,1);
	      else if (_orderedType == FERRARO_OUANGRAOUA)
		matching_distance = FerraroOuangraouaMatching(*_trees[i_tree],
							    *_trees[r_tree],
							    matching_sequence,_scale);
	      else if (_orderedType == FERRARO_OUANGRAOUA1)
		matching_distance = FerraroOuangraoua1Matching(*_trees[i_tree],
							    *_trees[r_tree],
							    matching_sequence,_scale);
	      else if (_orderedType == FERRARO_OUANGRAOUA2)
		matching_distance = FerraroOuangraouaMatching(*_trees[i_tree],
							    *_trees[r_tree],
							    matching_sequence,_scale);
	      else if (_orderedType == FERRARO_OUANGRAOUA2)
		matching_distance = FerraroOuangraouaMatching(*_trees[i_tree],
							    *_trees[r_tree],
							    matching_sequence,_scale);
	      else if (_orderedType == PARTIAL_ORDER)
		matching_distance = MatchPartialOrder(*_trees[i_tree],
							    *_trees[r_tree],
							    matching_sequence);
	    }
          if (_matchingType==TEST)
            {
              matching_distance=MatchTest(*_trees[i_tree],
                                          *_trees[r_tree],
                                          matching_sequence);
            }

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
          if ((_matchingType==BY_TOPOLOGY||_matchingType==ORDERED_MATCHING)&& _orderedType != FERRARO_OUANGRAOUA)
            {
              matching_sequence->putNbDel(tree_size1-seq_size);
              matching_sequence->putNbIns(tree_size2-seq_size);
            }
          //cout<<"MATCHING TIME "<<i_tree<<" with "<<r_tree<<" : "<<Time<<endl;
		  putTime(Time,i_tree,r_tree);
          //Time= 0.0;
        }
    }
}


//------------------------------------------------------------------------------------------------------
// LOCAL TOPOLOGICAL MATCHING
//------------------------------------------------------------------------------------------------------

void TreeMatch::localTopologicalMatching()
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

          Timer chrono;
          chrono.start();

          DistanceType matching_distance;
          if (_matchingType==BY_TOPOLOGY)
            {
              matching_distance=LocalMatchByTopology(*_trees[i_tree],
						     *_trees[r_tree],
						     matching_sequence);
            }
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
	  Sequence* s= getSequence(i_tree,r_tree);
	  cout<<"MATCHING TIME "<<i_tree<<" with "<<r_tree<<" : "<<Time<<endl;
          //Time= 0.0;
        }
    }
}


//------------------------------------------------------------------------------------------------------
// SEQUENCE MATCHING
//------------------------------------------------------------------------------------------------------
void TreeMatch::sequenceMatching()
{
  cout<<"Sequence matching"<<endl;
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

          Timer chrono;
          chrono.start();

          DistanceType matching_distance;
	      matching_distance=MatchSequences(*_trees[i_tree],
                                                *_trees[r_tree],
                                                matching_sequence);
          double Time= chrono.stop();

          putDistance(matching_distance,i_tree,r_tree);
          putSequence(matching_sequence,i_tree,r_tree);
//           int seq_size=matching_sequence->getSize();
//           int sub_number=0;
//           int mat_number=0;
//           matching_sequence->reset();
//           do
//             {
//               if (matching_sequence->getCurrent()->getCost()==0.0)
//                 {
//                   mat_number++;
//                 }
//               else
//                 {
//                   sub_number++;
//                 }
//             } while(matching_sequence->next());
//           matching_sequence->putNbMat(mat_number);
//           matching_sequence->putNbSub(sub_number);
// 	      matching_sequence->putNbDel(tree_size1-seq_size);
//           matching_sequence->putNbIns(tree_size2-seq_size);
          cout<<"MATCHING TIME "<<i_tree<<" with "<<r_tree<<" : "<<Time<<endl;
          //Time= 0.0;
        }
    }
}
//--------------------------------------------------------------------------------------
// DISTANCE FOR THE WAGNER AND FISHER'S ALGORITHM
//--------------------------------------------------------------------------------------
DistanceType TreeMatch::MatchSequences(TreeGraph& Tree1,
                                        TreeGraph& Tree2,
                                        Sequence * S)
{

  //  cout<<" Match by Topology "<<endl;
  // NodeCostFunction used
  NodeCost* MCF=new NodeCost(TOPOLOGIC);
  MatchingSequence* M=new MatchingSequence(Tree1,Tree2,*MCF,2);
  //effective matching : the result is stored in D
  DistanceType D=M->match();
  Sequence* s=new Sequence();
//   M->TreeList(0,0,*s);
//   S->putInsCost(M->getInsertCost());
//   S->putDelCost(M->getDeleteCost());
//   S->putSubCost(M->getSubstitutionCost());

//   if (s->getSize())
//     {
//       s->reset();
//       do
//         {
//           TreeNode* tree_node1=Tree1.getNode(s->getCurrent()->getIV());
//           TreeNode* tree_node2=Tree2.getNode(s->getCurrent()->getRV());
//           S->append(tree_node1->getVertex(),tree_node2->getVertex(),MCF->getChangingCost(tree_node1,tree_node2));
//         } while(s->next());
//     }
  delete (NodeCost*) MCF;
  delete (Sequence*) s;
  delete (MatchingSequence*) M;
  cout<<"End Match by Topology "<<endl;
  return(D);
}

//------------------------------------------------------------------------------------------------------
// SEQUENCE MATCHING
//------------------------------------------------------------------------------------------------------
void TreeMatch::selkowMatching()
{
  //  cout<<"Selkow's matching"<<endl;
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

          Timer chrono;
          chrono.start();

          DistanceType matching_distance;
	  matching_distance=MatchSelkow(*_trees[i_tree],
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
          cout<<"MATCHING TIME "<<i_tree<<" with "<<r_tree<<" : "<<Time<<endl;
          //Time= 0.0;
        }
    }
}
//--------------------------------------------------------------------------------------
// DISTANCE FOR THE WAGNER AND FISHER'S ALGORITHM
//--------------------------------------------------------------------------------------
DistanceType TreeMatch::MatchSelkow(TreeGraph& Tree1,
                                        TreeGraph& Tree2,
                                        Sequence * S)
{

  //  cout<<" Match by Topology "<<endl;
  // NodeCostFunction used
  NodeCost* MCF=new NodeCost(TOPOLOGIC);
  Selkow* M=new Selkow(Tree1,Tree2,*MCF);
  //effective matching : the result is stored in D
  DistanceType D=M->match();
  Sequence* s=new Sequence();
//   M->TreeList(0,0,*s);
//   Sequence* s1=new Sequence();
  M->ForestList(1,1,*s);
  S->putInsCost(M->getInsertCost());
  S->putDelCost(M->getDeleteCost());
  S->putSubCost(M->getSubstitutionCost());
  TreeNode* tree_node1=Tree1.getNode(0);
  TreeNode* tree_node2=Tree2.getNode(0);
  S->append(tree_node1->getVertex(),tree_node2->getVertex(),MCF->getChangingCost(tree_node1,tree_node2));

  if (s->getSize())
    {
      s->reset();
      do
        {
          TreeNode* tree_node1=Tree1.getNode(s->getCurrent()->getIV());
          TreeNode* tree_node2=Tree2.getNode(s->getCurrent()->getRV());
          S->append(tree_node1->getVertex(),tree_node2->getVertex(),MCF->getChangingCost(tree_node1,tree_node2));
        } while(s->next());
    }
  delete (NodeCost*) MCF;
  delete (Sequence*) s;
  delete (Selkow*) M;
  cout<<"End Selkow's Matching "<<endl;
  return(D);
}

