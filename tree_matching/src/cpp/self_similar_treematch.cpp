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
VPTOOLS_USING(Timer)
//------------------------------------------------------------------------------------------------------
// SELF SIMILAR TOPOLOGICAL MATCHING
//------------------------------------------------------------------------------------------------------

void TreeMatch::selfSimilarityTopologicalMatching()
{
  int nb_tree=_roots->entries();
  _trees.resize(nb_tree);
  _maxNbVertex=0;
  _maxOrder=0;
  for (int i=0;i<nb_tree;i++)
    {
      VId root=_roots->at(i);
      //cerr<<"root = "<<root<<endl;
      _trees[i]=new TreeGraph(*_mtg,root,TOPO,ANY);
      _maxOrder = I_MAX(_maxOrder,_trees[i]->getOrder());
      _maxNbVertex=I_MAX(_trees[i]->getNbVertex(),_maxNbVertex);
    }
  int nbTree;
  if (nb_tree == 1)
    nbTree = _trees[0]->getNbVertex();
  else 
    nbTree =  _trees[0]->getNbVertex() + _trees[1]->getNbVertex();
  _time.resize(1);
  _time[0].resize(1,0.0);
  _distances.resize(nbTree);
  _sequences.resize(nbTree);
  for (int tree=0;tree<nbTree;tree++)
    {
       _sequences[tree]=new SequenceVector(nbTree-tree,(Sequence*) NULL);
      _distances[tree].resize(nbTree-tree,0.0);
    }

  Timer chrono;
  chrono.start();
  if (nb_tree==1)
    DistanceType matching_distance=MatchByTopology(*_trees[0],
						   *_trees[0]);
  else
    DistanceType matching_distance=MatchByTopology(*_trees[0],
						   *_trees[1]);
  double Time = chrono.stop();
  cout<<"MATCHING TIME : "<<Time<<endl;
  Time= 0.0;
}

//------------------------------------------------------------------------------------------------------
// SELF SIMILAR COMPLEX MATCHING
//------------------------------------------------------------------------------------------------------

void TreeMatch::selfSimilarityComplexMatching()
{
  int nb_tree=_roots->entries();
  _trees.resize(nb_tree);
  _maxNbVertex=0;
  _maxOrder=0;
  for (int i=0;i<nb_tree;i++)
    {
      VId root=_roots->at(i);
      //cerr<<"root = "<<root<<endl;
      _trees[i]=new TreeGraph(*_mtg,root,TOPO,ANY);
      _maxOrder = I_MAX(_maxOrder,_trees[i]->getOrder());
      _maxNbVertex=I_MAX(_trees[i]->getNbVertex(),_maxNbVertex);
    }
  int nbTree;
  if (nb_tree == 1)
    nbTree = _trees[0]->getNbVertex();
  else 
    nbTree =  _trees[0]->getNbVertex() + _trees[1]->getNbVertex();
  _distances.resize(nbTree);
  _sequences.resize(nbTree);
  for (int tree=0;tree<nbTree;tree++)
    {
       _sequences[tree]=new SequenceVector(nbTree-tree,(Sequence*) NULL);
      _distances[tree].resize(nbTree-tree,0.0);
    }

  Timer chrono;
  chrono.start();
  if (nb_tree==1)
    DistanceType matching_distance=MatchByComplex(*_trees[0],
						   *_trees[0]);
  else
    DistanceType matching_distance=MatchByComplex(*_trees[0],
						   *_trees[1]);
  double Time = chrono.stop();
  cout<<"MATCHING TIME : "<<Time<<endl;
  //putTime(Time,0,1);
}

//--------------------------------------------------------------------------------------
//  d'un arbre dans un autre
//--------------------------------------------------------------------------------------
DistanceType TreeMatch::MatchEndSpaceFreeByTopology( TreeGraph& Tree1,
						     TreeGraph& Tree2)
{

  cout<<"Topological Begin-Space-Free Matching  (of one tree)"<<endl;
  // NodeCostFunction used
  NodeCost* MCF=new NodeCost(TOPOLOGIC);
  // making of the matching
  int optimization = 1;

  EndSpaceFreeMatching* M=new EndSpaceFreeMatching(Tree1,Tree2,*MCF );
  // effective matching : the result is stored in D
  DistanceType D=M->match_begin_T2();
  int size1 = Tree1.getNbVertex();
  int size2 = Tree2.getNbVertex();
  // making of the matching list
//  cerr << "\x0d" <<"Already Saved : 0% ...   " << flush;
  int nb_tree = _roots->entries();
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
  cout<<endl<<"Distance between Input Root = "<<M->getI_v()<<" and Reference Root "<<M->getR_v()<<" = "<<D<<endl;
  cerr<<"\x0d"<<flush;
  delete (NodeCost*) MCF;
  delete (EndSpaceFreeMatching*) M;
  return(D);
}


//------------------------------------------------------------------------------------------------------
// EndSpace Free TOPOLOGICAL MATCHING
//------------------------------------------------------------------------------------------------------


void TreeMatch::endSpaceFreeTopologicalMatching()
{
  int nb_tree=_roots->entries();
  _trees.resize(nb_tree);
  _maxNbVertex=0;
  _maxOrder=0;
  for (int i=0;i<nb_tree;i++)
    {
      VId root=_roots->at(i);
      //cerr<<"root = "<<root<<endl;
      _trees[i]=new TreeGraph(*_mtg,root,TOPO,ANY);
      _maxOrder = I_MAX(_maxOrder,_trees[i]->getOrder());
      _maxNbVertex=I_MAX(_trees[i]->getNbVertex(),_maxNbVertex);
    }
  cout<<"Entrer un numéro de vertex : ";
  int subtree;
  cin>>subtree;
//   _trees[0]->delSubTree(subtree);
//   cout<<"*****************"<<endl;
//   _trees[0]->print();
//   cout<<"*****************"<<endl;
  cout<<"Entrer un numéro de vertex : ";
  int subtree1;
  cin>>subtree1;
  _trees[1]->addSubTree(_trees[0],subtree,subtree1);
  cout<<"*****************"<<endl;
  _trees[1]->print();
  cout<<"*****************"<<endl;

  int nbTree;
  nbTree =  _trees[0]->getNbVertex() + _trees[1]->getNbVertex();
  _time.resize(1);
  _time[0].resize(1,0.0);
  _distances.resize(nbTree);
  _sequences.resize(nbTree);
  for (int tree=0;tree<nbTree;tree++)
    {
       _sequences[tree]=new SequenceVector(nbTree-tree,(Sequence*) NULL);
      _distances[tree].resize(nbTree-tree,0.0);
    }

  Timer chrono;
  chrono.start();
  DistanceType matching_distance=MatchEndSpaceFreeByTopology(*_trees[0],
							     *_trees[1]);
  double Time = chrono.stop();
  cout<<"MATCHING TIME : "<<Time<<endl;
  Time= 0.0;
}


//------------------------------------------------------------------------------------------------------
// FRACTAL TOPOLOGICAL MATCHING
//------------------------------------------------------------------------------------------------------
void TreeMatch::fractalTopologicalMatching()
{   
  int nb_tree=3;

  _trees.resize(nb_tree);
  _maxNbVertex=0;
  _maxOrder=0;
  
  VId root=_roots->at(0);
  //cerr<<"root = "<<root<<endl;
  _trees[0]=new TreeGraph(*_mtg,root,TOPO,ANY);
  //  _trees[0]->mtg_write(cout);
  _maxOrder = I_MAX(_maxOrder,_trees[0]->getOrder());
  _maxNbVertex=I_MAX(_trees[0]->getNbVertex(),_maxNbVertex);

  //cout<<"###########################################"<<endl;
  cout<<"Topological Fractal Matching "<<endl;

  // Récupération du numéro du vertex :
  int current_son_model = _trees[0]->getNumber(_roots->at(1));
  int current_root = _trees[0]->father(current_son_model);
  // cout<<"Current root = "<<current_root<<" - Current son = "<<current_son_model<<endl;
  // NodeCostFunction used
  NodeCost* MCF=new NodeCost(TOPOLOGIC);
  // making of the matching
  
  int nbTree = 1;
  _distances.resize(nb_tree);
  _sequences.resize(nb_tree);
  for (int tree=0;tree<nb_tree;tree++)
    {
       _sequences[tree]=new SequenceVector(nb_tree-tree,(Sequence*) NULL);
      _distances[tree].resize(nb_tree-tree,0.0);
    }

  Timer chrono;
  chrono.start();
  
  DistanceType distance_to_self_similar_tree = 0.0;
  
  while (current_son_model !=0){
      
    //list<int>::const_iterator begin_children =  (_trees[0]->sons(current_root))->begin();
    // On récupère le premier fils
    
    _trees[1] = new TreeGraph(_trees[0],current_son_model);

    //cout<<"###### Tree 1 ###################"<<endl;
    //_trees[1]->print();
    int i = 1;
    int i_child = 1;
    int nb_sons=(_trees[0]->sons(current_root))->size();
    //cout<<"nombre de fils = "<<nb_sons<<endl;
    while (i<=nb_sons){
      int current_child = (_trees[0]->child(current_root,i_child));
      //cout<<"nombre de fils = "<<(_trees[0]->sons(current_root))->size()<<endl;
      //cout<<i<<"eme son : "<<current_child<<endl;
      i++;
      i_child++;
      //cout<<"Father = "<<current_root<<" - Son Model = "<<current_son_model<<" - Current son =  "<<current_child<<endl;
      //      if (current_child != current_son_model){
	_trees[2] = new TreeGraph(_trees[0],(current_child));
	
	//cout<<"###### Tree 2 ###################"<<endl;
	//_trees[2]->print();
	EndSpaceFreeMatching* M=new EndSpaceFreeMatching(*_trees[1],*_trees[2],*MCF );

	// effective matching : the result is stored in D
	DistanceType D=M->match_begin_T2();

	cout<<endl<<"Distance between Input Root = "<<(_trees[1]->getTreeNode(M->getI_v()))->getVertex()
	   <<" and Reference Root "<<(_trees[2]->getTreeNode(M->getR_v()))->getVertex()<<" = "<<D<<endl;
	if (D>0){
	  _trees[0]->delSubTree(current_child);
	  //i_child-- ;
	  //_trees[0]->print();
	  //cout<<"###############################"<<endl;
	  //cout<<current_root<<" - "<<M->getI_v()<<endl;

	  _trees[0]->addSubTree(_trees[1],M->getI_v(),current_root);
	  _trees[0]->print();
	}
	distance_to_self_similar_tree += D;
	delete (EndSpaceFreeMatching*) M;
	// }
	//      begin_children++;
    }
    current_son_model = current_root;
    current_root = _trees[0]->father(current_son_model);
  }
  putDistance(distance_to_self_similar_tree,0,1);
  cout<<"Self Similar value = " << distance_to_self_similar_tree <<endl;
  double Time = chrono.stop();
  cout<<"MATCHING TIME : "<<Time<<endl;
  Time= 0.0;
  delete (NodeCost*) MCF;
}

