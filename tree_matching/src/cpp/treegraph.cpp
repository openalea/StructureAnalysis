/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2000 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): P.ferraro (pascal.ferraro@cirad.fr)
 *
 *       $Source: /usr/cvsmaster/AMAPmod/src/TreeMatching/treegraph.cpp,v $
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


#include "treegraph.h"

using namespace std;

void TreeGraph::makeNode(VId vertex,int father,VId compo_father,TreeType type)
{
  TreeNode* tnode=new TreeNode(*_mtg);
  NodeList* nlist=new NodeList();

  VIdList* vertex_sons;
  VId plant_root;
  switch(type)
    {
    case TOPO   :{
      vertex_sons=_mtg->topoSons(vertex,_edge);
      plant_root = _mtg->findGlobalTopoRoot(vertex, ANYTOPO);
    }break;
    case COMPO  : {
      vertex_sons=_mtg->localTopoSons(vertex,_edge);  
      plant_root = _mtg->findLocalTopoRoot(vertex, ANYTOPO);
    }break;
    default     : assert(0);break;
    }


  //  VId plant_root = _mtg->findGlobalTopoRoot(vertex, ANYTOPO);
  int order = _mtg->pathEdgeNb(plant_root,vertex,PLUS);
  _order = I_MAX(_order,order);
  tnode->putOrder(order);


  int nb_child=vertex_sons->entries();

  if (nb_child>_degree) _degree=nb_child;

  if (father!=-1)
    {
      putNode(father,_number);
//       int depth=_treenodes.at(father)->depth()+1;
//       list<TreeNode*>::iterator begin;
//       begin = _treenodes.begin();
//       for (int i=0;i<father;i++)
//         begin++;
//       int depth=(*begin)->depth()+1;
      int depth=(_treenodes[father])->depth()+1;
      tnode->putDepth(depth);
      _depth=I_MAX(_depth,depth);
    }
  else
    {
      tnode->putDepth(1);
      _depth=1;
    }
  
  tnode->putFather(father);
  tnode->putNumber(_number);
  int num = _number;
  tnode->putVertex(vertex);
  tnode->putValue(1);
  tnode->putComplex(_mtg->compoFather(vertex));

  
  // On entre les valeurs correspondant aux attributs
  if (_valued)
    {
      tnode->resize(_functions->size());
      for (int i=0;i<_functions->size();i++)
        {
          NodeFunctionList::iterator begin;
          begin = _functions->begin();
          for (int j=0;j<i;j++)
            begin++;
          tnode->putValue(i,(*begin)(vertex));
        }
    }
  putTreeNode(tnode);
  putNodeList(nlist);
  _number++;

  VIdListIter next(*vertex_sons);
  VId son;
  int numSon;
  if(next()){ 
    son = next.key ();
    NodeList* nlist1=new NodeList();
    putClassNodeList(nlist1 );
    numSon = _number;
    putNodeinClass(_numClass++,numSon);
    _class[numSon] = _numClass-1;
    makeNode(son,tnode->getNumber(),compo_father,type);
  }
  
  while(next())
    {
      VId prec_son = son;
      int prec_numSon = numSon; 
      son = next.key();
      numSon = _number;
      if(_mtg->existsFName("ind")){
	if(_mtg->si_feature(son,"ind")->i == _mtg->si_feature(prec_son,"ind")->i)
	  {
	    putNodeinClass(_class[prec_numSon],numSon);
	    _class[numSon] = _class[prec_numSon];
	  }
	else
	  {
	    NodeList* nlist1=new NodeList();
	    putClassNodeList(nlist1);
	    putNodeinClass(_numClass++,numSon);
	    _class[numSon] = _numClass-1;
	  }
      }else{
	putNodeinClass(_class[prec_numSon],numSon);
	_class[numSon] = _class[prec_numSon];
      }
      
      makeNode(son,tnode->getNumber(),compo_father,type);
      
    }
  delete (VIdList*) vertex_sons;
  
  tnode->putNumPostfix(_numPostfix);
  if(_postfix.size() < _numPostfix+1)
    _postfix.resize(_numPostfix+1);
  _postfix[_numPostfix]= num;
  _numPostfix++;
  
}

void TreeGraph::makeNode(VId vertex,int father,VId compo_father,TreeType type,VIdList* vlist)
{


 TreeNode* tnode=new TreeNode(*_mtg);
  NodeList* nlist=new NodeList();

  VIdList* vertex_sons;
  VId plant_root;
  switch(type)
    {
    case TOPO   :{
      vertex_sons=_mtg->topoSons(vertex,_edge);
      plant_root = _mtg->findGlobalTopoRoot(vertex, ANYTOPO);
    }break;
    case COMPO  : {
      vertex_sons=_mtg->localTopoSons(vertex,_edge);  
      plant_root = _mtg->findLocalTopoRoot(vertex, ANYTOPO);
    }break;
    default     : assert(0);break;
    }


  //  VId plant_root = _mtg->findGlobalTopoRoot(vertex, ANYTOPO);
  int order = _mtg->pathEdgeNb(plant_root,vertex,PLUS);
  _order = I_MAX(_order,order);
  tnode->putOrder(order);


  
  if (father!=-1)
    {
      putNode(father,_number);
      cerr<<_treenodes.size()<<" => father ="<<father<<endl;
      int depth=(_treenodes[father])->depth()+1;
//       list<TreeNode*>::iterator begin;
//       begin = _treenodes.begin();
//       for (int i=0;i<father;i++)
//         begin++;
//       int depth=(*begin)->depth()+1;
      tnode->putDepth(depth);
      _depth=I_MAX(_depth,depth);
    }
  else
    {
      tnode->putDepth(1);
      _depth=1;
    }
  
  tnode->putFather(father);
  tnode->putNumber(_number);
  int num = _number;
  tnode->putVertex(vertex);
  tnode->putValue(1);
  tnode->putComplex(_mtg->compoFather(vertex));

  
  // On entre les valeurs correspondant aux attributs
  if (_valued)
    {
      tnode->resize(_functions->size());
      for (int i=0;i<_functions->size();i++)
        {
          NodeFunctionList::iterator begin;
          begin = _functions->begin();
          for (int j=0;j<i;j++)
            begin++;
          tnode->putValue(i,(*begin)(vertex));
        }
    }
  putTreeNode(tnode);
  putNodeList(nlist);
  _number++;
  

  int nb_child=0;


  VIdListIter next(*vertex_sons);
  VId son=-1;
  int numSon;
  int first=0;
  int b=0;
  while(next()){
	VId prec_son = son;
	int prec_numSon = numSon; 
	son = next.key();

	if(first==0){	
      b=0;
	  VIdListIter next2(*vlist);
	  VId v;    
	  while(next2()){ 
		v = next2.key();
		if((_mtg->isTopoDescendantOf(v,son,ANY))&&(b==0)){b=1;first=1;}
	  }
    
    
	  if(b==1){
		nb_child++;
		NodeList* nlist1=new NodeList();
        putClassNodeList(nlist1);
        numSon = _number;
        putNodeinClass(_numClass++,numSon);
        _class[numSon] = _numClass-1;
        makeNode(son,tnode->getNumber(),compo_father,type,vlist);
	  }
	} 
    else
    {  
	  b=0;
	  VIdListIter next2(*vlist);
      VId v;    
      while(next2()){ 
		v = next2.key();
		if((_mtg->isTopoDescendantOf(v,son,ANY))&&(b==0)){b=1;}
	  }
      

      if(b==1){
		nb_child++;
	    numSon = _number;
	    if(_mtg->existsFName("ind")){
	      if((_mtg->si_feature(son,"ind"))->i == (_mtg->si_feature(prec_son,"ind"))->i){
	        putNodeinClass(_class[prec_numSon],numSon);
	        _class[numSon] = _class[prec_numSon];
		  }
	      else
		  {
	        NodeList* nlist1=new NodeList();
	        putClassNodeList(nlist1);
	        putNodeinClass(_numClass++,numSon);
	        _class[numSon] = _numClass-1;
		  }
		}
		else
		{
	      putNodeinClass(_class[prec_numSon],numSon);
	      _class[numSon] = _class[prec_numSon];
		}
	    makeNode(son,tnode->getNumber(),compo_father,type,vlist);
	  }  
	}
  }
  
  if (nb_child>_degree) _degree=nb_child;

  delete (VIdList*) vertex_sons;
  
  tnode->putNumPostfix(_numPostfix);
  if(_postfix.size() < _numPostfix+1)
    _postfix.resize(_numPostfix+1);
  _postfix[_numPostfix]= num;
  _numPostfix++;
 
  }

void TreeGraph::makeComplexNode(VId vertex,int father,VId compo_father,TreeType type)
{
  TreeNode* tnode=new TreeNode(*_mtg);
  NodeList* nlist=new NodeList();
  
  VIdList* vertex_sons;
  switch(type)
    {
    case TOPO 	: vertex_sons=_mtg->topoSons(vertex,_edge);break;
    case COMPO 	: vertex_sons=_mtg->localTopoSons(vertex,_edge);break;
    default 	: assert(0);break;
    }
  
  
  VId plant_root = _mtg->findGlobalTopoRoot(vertex, ANYTOPO);
  int order = _mtg->pathEdgeNb(plant_root,vertex,PLUS);
  _order = I_MAX(_order,order);
  tnode->putOrder(order);
  
  
  int nb_child=vertex_sons->entries();
  
  if (nb_child>_degree) _degree=nb_child;
  
  if (father!=-1) 
    {
      putNode(father,_number);      
      int depth=_treenodes[father] ->depth()+1;
//       list<TreeNode*>::iterator begin;
//       begin = _treenodes.begin();
//       for (int i=0;i<father;i++)
// 	begin++;
//       int depth=(*begin)->depth()+1;
      tnode->putDepth(depth);
      _depth=I_MAX(_depth,depth);
    }
  else
    {
      tnode->putDepth(1);
      _depth=1;	
    }
  
  tnode->putFather(father);
  tnode->putNumber(_number);
  int num = _number;
  tnode->putVertex(vertex);
  tnode->putValue(1);
  tnode->putComplex(_mtg->compoFather(vertex));
  
  // On entre les valeurs correspondant aux attributs
  if (_valued)
    {
      tnode->resize(_functions->size());
      for (int i=0;i<_functions->size();i++)
	{
	  NodeFunctionList::iterator begin;
	  begin = _functions->begin();
	  for (int j=0;j<i;j++)
	    begin++;
 	  tnode->putValue(i,(*begin)(vertex));
	}
    }
  putTreeNode(tnode);
  putNodeList(nlist);
  _number++;
  
  VIdListIter next(*vertex_sons);
  while(next())
    {
      VId son = next.key();
      if (fatherIsInComplex(_number-1))
	makeComplexNode(son,tnode->getNumber(),compo_father,type);
    }
  delete (VIdList*) vertex_sons;
  tnode->putNumPostfix(_numPostfix);
  if(_postfix.size() < _numPostfix+1)
    _postfix.resize(_numPostfix+1);
  _postfix[_numPostfix]= num;
  _numPostfix++;
  
}

int TreeGraph::getNumber(int vertex) const
{
//   vector<TreeNode*>::const_iterator begin;
//   begin = _treenodes.begin();
//   int i = 0;
//   while (((*begin)->getVertex()!=vertex)&&(i<_nbVertex)){
//     begin++;
//     i++;
//   }
  if (vertex>=_nbVertex)
    return -1;
  else
    return((_treenodes[vertex])->getNumber());
}

void TreeGraph::makeAxeNode(VId vertex,int father,int MAX_ORDER,int position)
{
  TreeNode* axe_node=new TreeNode(*_mtg);
  NodeList* nlist=new NodeList();

  // each treenode represents an axis

  VId plant_root = _mtg->findGlobalTopoRoot(vertex, ANYTOPO);
  int order = _mtg->pathEdgeNb(plant_root,vertex,PLUS);
  _order = I_MAX(_order,order);
  axe_node->putOrder(order);

  if (father!=-1)
    {
      putNode(father,_number);
    }

  int number=_number;
  axe_node->putPosition(position);
  axe_node->putOrder(order);
  axe_node->putFather(father);
  axe_node->putNumber(number);
  axe_node->putVertex(vertex);
  axe_node->putComplex(_mtg->compoFather(vertex));
  axe_node->putValue(1);
  if (father!=-1)
    {
      // On identifie l'élément n° father de _treenode
//       list<TreeNode*>::iterator begin;
//       begin = _treenodes.begin();
//       for (int i=0;i<father;i++)
//         begin++;
//       int depth=(*begin)->depth()+1;
      int depth=(_treenodes[father])->depth()+1;

      axe_node->putDepth(depth);
      if (_depth<depth) _depth=depth;
    }
  else
    {
      axe_node->putDepth(1);
      _depth=1;
    }
  if (_valued)
    {
      axe_node->resize(_functions->size());
      for (size_t i=0;i<_functions->size();i++)
        {
          NodeFunctionList::iterator begin;
          begin = _functions->begin();
          for (size_t j=0;j<i;j++)
            begin++;
          axe_node->putValue(i,(*begin)(vertex));
          // axe_node->putValue(i,(*_functions->at(i))(vertex));
        }
    }
  putTreeNode(axe_node);
  putNodeList(nlist);

  VIdList* axe=_mtg->findLessSequence(vertex);
  _nbAxisVertex=_nbAxisVertex+axe->entries();
  // We explore the axis as long as their order is less than MAX_ORDER
  if (order<MAX_ORDER)
    {
      // getting of the axis
      VIdListIter next_on_axe(*axe);
      int position_on_axe=1;
      while(next_on_axe())
        {
          VId axe_vertex=next_on_axe.key();
          // Getting of the vertices corresponding to new axises
          VIdList* sons=_mtg->topoSons(axe_vertex,PLUS);
          VIdListIter next_son(*sons);
          while(next_son())
            {
              _number++;
              makeAxeNode(next_son.key(),number,MAX_ORDER,position_on_axe);
            }
          delete (VIdList*) sons;
          position_on_axe++;
        }
    }
  delete (VIdList*) axe;
  int nb_child=nlist->size();
  if (nb_child>_degree) _degree=nb_child;

}

TreeGraph::TreeGraph( )
{
  _mtg=0;
  _degree=0;
  _depth=0;
  _number=0;
 _numPostfix=0;
  _root=0;
  _order=0;
  _nbVertex=0;
  _nbAxisVertex=0;
  _valued=FALSE;
  _functions = (NodeFunctionList*) NULL;
  _mtg = (MTG*) NULL;

}


TreeGraph::TreeGraph(MTG& mtg,VId root,TreeType type,EType edge_type,AmlBoolean valued,NodeFunctionList* functions)
{
  // MAKING OF THE TREE
  _mtg=&mtg;
  _degree=1;
  _depth=0;
  _number=0;
 _numPostfix=0;
  _numClass=0;
  int father=-1;
  _order=0;
  _edge=edge_type;
  _valued=valued;
  _functions=functions;
  switch (type)
    {
    case TOPO   :_root=root;break;
    case COMPO  :_root=mtg.findFirstComponentAtScale(root,mtg.vscale(root)+1);break;
    default             :assert(0);break;
    }
  _class.resize(_mtg->vertexNb()+2);
  NodeList* nlist1=new NodeList();
  putClassNodeList(nlist1 );
  putNodeinClass(_numClass,0);
  _class[0]=_numClass++;
  makeNode(_root,father,root,type);
  //    if (_valued) normalizeNodes();
  _nbVertex=_number;
_nbClass=_numClass;
}

int TreeGraph::getChildLess(int vertex) const
{
  int order = getTreeNode(vertex)->getOrder();
  int child_number=getNbChild(vertex);
//   list<NodeList*>::const_iterator begin;
//   begin = _outNodeTable.begin();
//   for (int i=0;i<vertex;i++)
//     begin++;
   vector<int>::const_iterator begin_child;
   begin_child = (_outNodeTable[vertex])->begin();
  int j= 0;
  for ( j=0;(getTreeNode(*begin_child)->getOrder()!=order)&&j<child_number;j++)
    begin_child++;
  if (j==child_number)
    return -1;
  else
    return(*begin_child);
}

void TreeGraph::printNodedClass() const
{

  vector<NodeList*>::const_iterator begin;
  begin = _classNodeTable.begin();
  for (int i=0;i<getNbClass();i++){
    cout<<"Class "<<i<<" :";
    vector<int>::const_iterator begin_node;
    begin_node = (*begin)->begin();
    while(begin_node !=  (*begin)->end()){
      cout<<(*begin_node)<<" ; ";
      begin_node++;
    }
    cout<<endl;
    begin++;
  }
}


/* Making of tree Graphs on a given complex */

TreeGraph::TreeGraph(MTG& mtg,VId root,TreeType type,AmlBoolean valued,NodeFunctionList* functions)
{
  // MAKING OF THE TREE	
  _mtg=&mtg; 
  _degree=1;
  _depth=0;
  _number=0;
_numPostfix=0;
   int father=-1;
  _order=0;
  _edge=ANY;
  _valued=valued;
  _functions=functions;
  switch (type)
    {
    case TOPO 	:_root=root;break;
    case COMPO 	:_root=mtg.findFirstComponentAtScale(root,mtg.vscale(root)+1);break;
    default		:assert(0);break;
    }
  makeComplexNode(_root,father,root,type);
  //	if (_valued) normalizeNodes();
  _nbVertex=_number;
}

TreeGraph::TreeGraph(MTG& mtg,VIdList* vlist, TreeType type ,EType edge_type,AmlBoolean valued,NodeFunctionList* functions){
  // MAKING OF THE TREE
  _mtg=&mtg;
  _degree=1;
  _depth=0;
  _number=0;
  _numPostfix=0;
  _numClass=0;
  int father=-1;
  _order=0;
  _edge=edge_type;
  _valued=valued;
  _functions=functions;
  
  //SEARCH THE ROOT


  VId root;
  VId v;
  VIdListIter next(*vlist);
  if(next()){root=next.key ();}
 
  while(next()){ 
    v = next.key ();
    if(v<root){root=v;}
  }
 
  
  switch (type)
    {
    case TOPO   :_root=root;break;
    case COMPO  :_root=mtg.findFirstComponentAtScale(root,mtg.vscale(root)+1);break;
    default             :assert(0);break;
    }
  _class.resize(_mtg->vertexNb()+2);
  NodeList* nlist1=new NodeList();
  putClassNodeList(nlist1 );
  putNodeinClass(_numClass,0);
  _class[0]=_numClass++;
  makeNode(_root,father,root,type,vlist);
  //    if (_valued) normalizeNodes();
  _nbVertex=_number;
  _nbClass=_numClass;
}


TreeGraph::TreeGraph(MTG& mtg,VId root,int MAX_ORDER,AmlBoolean valued,NodeFunctionList* functions)
{
  // MAKING OF THE TREE
  _mtg=&mtg;
  _degree=1;
  _depth=0;
  _number=0;
_numPostfix=0;

  _root=root;
  _order=0;
  int position=0;
  int father=-1;
  _nbAxisVertex=0;
  _edge=ANY;
  _valued=valued;
  _functions=functions;
  makeAxeNode(root,father,MAX_ORDER,position);
  //    if (_valued) normalizeNodes();
  _nbVertex=_number+1;
}

TreeGraph::TreeGraph(char* fich_name, MTG& mtg, VId root, AmlBoolean valued,NodeFunctionList* functions)
{
  string line;
  ifstream inFile(fich_name,ios::in);
  int NbChild=0;

  //line.readLine(inFile,TRUE);

  size_t lineSize = 1024;
  char * lineBuf = new char[lineSize];
  inFile.getline(lineBuf,lineSize);
  line = string(lineBuf);

  _nbVertex=atoi(line.data());
  _degree=1;
  _depth=1;
  _edge=ANY;
  _root=root;
  _valued=valued;
  _mtg=&mtg;
  _order=0;
  _functions=functions;
  _number=_nbVertex;
_numPostfix=0;
// FOR EACH NODE
  for(int i=0;i<_nbVertex;i++)
    {

      // READING OF THE NUMBER OF CHILD
      inFile.getline(lineBuf,lineSize);
      line = string(lineBuf);

      NbChild=atoi(line.data());

      // COMPUTING OF THE DEGREE
      if (NbChild>_degree) { _degree=NbChild; };

      // MAKING OF THE CHILD LIST
      TreeNode* node=new TreeNode(*_mtg);
      NodeList* nlist=new NodeList();

      // READING OF THE NODE VALUE
      inFile.getline(lineBuf,lineSize);
      line = string(lineBuf);
      node->putDepth(atoi(line.data()));
      _depth=I_MAX(_depth,atoi(line.data()));

      inFile.getline(lineBuf,lineSize);
      line = string(lineBuf);
      node->putFather(atoi(line.data()));

      inFile.getline(lineBuf,lineSize);
      line = string(lineBuf);
      node->putNumber(atoi(line.data()));

      inFile.getline(lineBuf,lineSize);
      line = string(lineBuf);
      node->putValue(atoi(line.data()));

      if (_valued)
        {
          node->resize(_functions->size());
          for (size_t j=0;j<_functions->size();j++)
            {
              NodeFunctionList::iterator begin;
              begin = _functions->begin();
              for (size_t k=0;k<j;k++)
                begin++;
              node->putValue(j,(*begin)(root));
            }
        }

      node->putVertex(root);
      node->putComplex(root);

      putTreeNode(node);
      putNodeList(nlist);

      // READING OF CHILDS' NUMBER
      for(int k=1;k<=NbChild;k++)
        {
          inFile.getline(lineBuf,lineSize);
          line = string(lineBuf);
          putNode(i,atoi(line.data()));
        }
    }
  //    if (_valued) normalizeNodes();
}


TreeGraph::TreeGraph(TreeGraph* T,int vertex)
{

    // MAKING OF THE TREE
  _mtg=T->getMTG();
  _degree=1;
  int root_depth = (T->getTreeNode(vertex))->depth();
  _depth= T->getDepth()-root_depth+1;
  _number=0;
  _numPostfix=0;

  _root=(T->getTreeNode(vertex))->getVertex();
  _order=0;
  int position=0;
  int father=-1;
  _edge=ANY;
  // Not implemented
//   _valued=valued;
//   _functions=functions;


  vector<TreeNode*>::iterator itreenode;
  vector<NodeList*>::iterator i_node_table;
  itreenode = _treenodes.begin();
  i_node_table = _outNodeTable.begin();
  int i;
  int nb_desc = T->getNbDesc(vertex)+1;
  _nbVertex=nb_desc;
  for (i=vertex;i<vertex+nb_desc;i++){
    TreeNode* tnode = new TreeNode(*(T->getTreeNode(i)));
    tnode->putNumber(i-vertex);
    tnode->putDepth(tnode->depth()-root_depth+1);
    if (i-vertex>0)
      tnode->putFather(tnode->father()-vertex);
    else 
      tnode->putFather(-1);
    putTreeNode(tnode);
    // On cosntruit la liste des fils
    NodeList* new_sonlist = new NodeList();
    NodeList sonlist =*(T->sons(i));
    vector<int>::iterator it_son_list;
    it_son_list = sonlist.begin();
    while (it_son_list != sonlist.end()){
      new_sonlist->push_back((*it_son_list) - vertex);
      it_son_list++;
    }
    if (_degree<new_sonlist->size())
      _degree = new_sonlist->size();
    putNodeList(new_sonlist);    
  }
}




//Destructor
TreeGraph::~TreeGraph()
{
  for (vector<NodeList*>::iterator i = _outNodeTable.begin(); i != _outNodeTable.end(); ++i)
    delete *i;

  for (vector<TreeNode*>::iterator i = _treenodes.begin(); i != _treenodes.end(); ++i)
    delete *i;
  
  _treenodes.clear();
  _outNodeTable.clear();
}

int TreeGraph::isLeaf(int node) const
{
  return((getNbChild(node)==0));
}


int TreeGraph::father(int node) const
{
//   list<TreeNode*>::const_iterator begin;
//   begin = _treenodes.begin();
//   for (int i=0;i<node;i++)
//     begin++;
//   return((*begin)->father());
  return (_treenodes[node])->father();
}

int TreeGraph::child(int node,int child_number) const
{
  if (child_number>getNbChild(node))
    {
      return(-1);
    }
  else
    {
      return ((*_outNodeTable[node])[child_number-1]);
//       list<NodeList*>::const_iterator begin;
//       begin = _outNodeTable.begin();
//       for (int i=0;i<node;i++)
//         begin++;
//       list<int>::const_iterator begin_child;
//       begin_child = (*begin)->begin();
//       for (int j=0;j<child_number-1;j++)
//         begin_child++;
//       return(*begin_child);
      //      return((*begin)->at(child_number-1));
    }
}


 int TreeGraph::getRealNumber(int postfix) const { return(_postfix[postfix]);}

int  TreeGraph::leftBrother(int node) const
{
  int fat = father(node);
//   list<NodeList*>::const_iterator begin;
//   begin = _outNodeTable.begin();
//   for (int i=0;i<fat;i++)
//     begin++;
  vector<int>::const_iterator begin_brother;
  begin_brother = (_outNodeTable[fat])->begin();
  int left = -1;
  while ((*begin_brother)!=node){
    left = *begin_brother;
    begin_brother++;
  }
  return (left) ;
}

int  TreeGraph::rightBrother(int node) const
{
  int fat = father(node);
//   int fat = father(node);
//   list<NodeList*>::const_iterator begin;
//   begin = _outNodeTable.begin();
//   // We get the set of sons of fat
//   for (int i=0;i<fat;i++)
//     begin++;
  vector<int>::const_iterator begin_brother;
  begin_brother = (_outNodeTable[fat])->end();
  int right = -1;
  begin_brother--;
  while ((*begin_brother) != node){
    right = *begin_brother;  
    begin_brother--;
  }
    begin_brother++;
    right = *begin_brother;  
  return (right) ;
}


const NodeList* TreeGraph::sons(int vertex) const
{
//   list<NodeList*>::const_iterator begin;
//   begin = _outNodeTable.begin();
//   for (int i=0;i<vertex;i++)
//     begin++;
//   return(*begin);
  return(_outNodeTable[vertex]);
}

void TreeGraph::removeSon(int vertex,int son,int decalage)
{
  vector<int>::iterator i_node_list;
  i_node_list = (_outNodeTable[vertex])->begin();
  while ((i_node_list != (_outNodeTable[vertex])->end()) && ((*i_node_list) != son ))
    i_node_list++;
  vector<int>::iterator n_node_list=i_node_list;
   while (n_node_list != (_outNodeTable[vertex])->end()){
     (*n_node_list) -= decalage;
     n_node_list++;
   }
  (_outNodeTable[vertex])->erase(i_node_list);
}

void TreeGraph::shiftSonNumber(int vertex,int son,int decalage)
{
  vector<int>::iterator i_node_list;
  i_node_list = (_outNodeTable[vertex])->begin();
  while ((i_node_list != (_outNodeTable[vertex])->end()) && ((*i_node_list) != son ))
    i_node_list++;
  i_node_list++;
  vector<int>::iterator n_node_list=i_node_list;
  while (n_node_list != (_outNodeTable[vertex])->end()){
    (*n_node_list) += decalage;
    n_node_list++;
  }
}


TreeNode* TreeGraph::getTreeNode(int vertex) const
{
//   list<TreeNode*>::const_iterator begin;
//   begin = _treenodes.begin();
//   for (int i=0;i<vertex;i++)
//     begin++;
//   return(*begin);
  return _treenodes[vertex];
}

int  TreeGraph::childIsInComplex(int vertex,int child_number) const
{
  if (child_number>getNbChild(vertex))
    {
      return(-1);
    }
  else
    {
      if (getTreeNode(vertex)->getComplex()== getTreeNode(child(vertex,child_number))->getComplex())
        //if (_treenodes.at(vertex)->getComplex()==(_treenodes.at(child(vertex,child_number)))->getComplex())
        return(1);
      else
        return (0);
    }
}

int  TreeGraph::childIsInAxis(int vertex,int child_number) const
{
  if (child_number>getNbChild(vertex))
    {
      return(-1);
    }
  else
    {
      if (getTreeNode(vertex)->getOrder()== getTreeNode(child(vertex,child_number))->getOrder())
        return(1);
      else
        return (0);
    }
}


int  TreeGraph::fatherIsInComplex(int vertex) const
{
  if (getTreeNode(vertex)->getComplex()==getTreeNode(father(vertex))->getComplex())
    return(1);
  else
    return (0);
}


int TreeGraph::getNbChild(int vertex) const
{
  if (vertex==-1)
    return(0);

  else
    return(sons(vertex)->size());
}

int TreeGraph::getNbDesc(int vertex) const
{
  if (isLeaf(vertex))
    return(0);
  else
    {
      int nb_desc= getNbChild(vertex);
      for (int i = 1; i<=getNbChild(vertex);i++)
        {
          nb_desc = nb_desc+getNbDesc(child(vertex,i));
        }
      return(nb_desc);
    }
}
int TreeGraph::getRoot() const
{
  return(_root);
}

int TreeGraph::getDegree() const
{
  return(_degree);
}

int TreeGraph::getOrder() const
{
  return(_order);
}

int TreeGraph::getDepth() const
{
  return(_depth);
}

int TreeGraph::getNbVertex() const
{
  return(_nbVertex);
}


int TreeGraph::getNbClass() const
{
  return(_nbClass);
}
int TreeGraph::getNbAxisVertex() const
{
  return(_nbAxisVertex);
}

TreeNode* TreeGraph::getNode(int i) const
{
  return(getTreeNode(i));
}

void TreeGraph::putNode(int node,int child)
{
  (_outNodeTable[node])->push_back(child);
}


void TreeGraph::putNodeinClass(int clas,int node)
{
  (_classNodeTable[clas])->push_back(node);
}

void TreeGraph::putTreeNode(TreeNode* tree_node)
{
  _treenodes.push_back(tree_node);
}

void TreeGraph::putNodeList(NodeList* node_list)
{
  _outNodeTable.push_back(node_list);
}

void TreeGraph::putClassNodeList(NodeList* node_list)
{
  _classNodeTable.push_back(node_list);
}
int TreeGraph::isNull()
{
  return ((_mtg==0)&&(_nbVertex==0)&&(_degree==0)&&(_depth==0));
}

void TreeGraph::print()
{
  cout<<" ROOT : "<<_root<< std::endl;
  cout<<" DEGREE : "<<_degree<<std::endl;
  cout<<" NB_VERTEX : "<<_nbVertex<<std::endl;
  char TAB='\t';
  for(int i=0;i<getNbVertex();i++)
    {
      // READING OF THE NODE VALUE
      TreeNode* node=getNode(i);
      node->print();
      cout<<TAB<<" CHILD : "<<endl;
      int n=0;
      for(int k=1;k<=getNbChild(i);k++)
        {
          cout<<TAB<<child(i,k);
          if (n!=10) {n++;} else {cout<<endl;n=0;};
        }
      cout<<endl;
    }
}

bool TreeGraph::mtg_write( const char *path ) const
{
  bool status = false;
  ofstream out_file(path);
  if (!out_file) {
    status = false;
    cerr<<"Can't open file"<<endl;
  }
  else {
    status = true;
    mtg_write(out_file);
  }
  return status;
}


ostream& TreeGraph::mtg_write(ostream &os) const
{
  vector<TreeNode*>::const_iterator itreenode;
  itreenode = _treenodes.begin();
  os<<"CODE:	FORM-A			"<<endl;
  os<<endl;
  os<<endl;
  os<<"CLASSES:				"<<endl;
  os<<"SYMBOL\tSCALE\tDECOMPOSITION\tINDEXATION\tDEFINITION"<<endl;
  os<<"$ \t0\tFREE\tFREE\tIMPLICIT"<<endl;
  os<<"P \t 1 \t CONNECTED \t FREE \t EXPLICIT"<<endl;
  os<<"E \t 2 \t <-LINEAR \t FREE \t EXPLICIT"<<endl;
  os<<endl;
  os<<"DESCRIPTION :				"<<endl;
  os<<"LEFT \t RIGHT \t RELTYPE \t MAX	"<<endl;
  os<<"E \t E \t < \t 1 	"<<endl;
  os<<"E \t E \t + \t ? 	"<<endl;
  os<<endl;
  os<<"FEATURES:   "<<endl;
  os<<"NAME \t TYPE"<<endl;
  os<<"VId \t INT"<<endl;
  os<<endl;
  os<<"MTG:"<<endl;
  os<<"ENTITY-CODE";
  int i;
  for (i=0;i<_depth;i++)
    os<<"\t";
  os<<"VId"<<endl;
  os<<"P1/";
  os<<"E"<<((*itreenode)->getNumber());
  for (i=0;i<_depth-(*itreenode)->depth()+1;i++)
    os<<"\t";
  os<<(*itreenode)->getVertex()<<endl;
  itreenode++;
  while (itreenode != _treenodes.end()){
    for (i=0;i<(*itreenode)->depth()-1;i++)
      os<<"\t";
    os<<"+E"<<(*itreenode)->getNumber();
    for (i=0;i<_depth-(*itreenode)->depth()+1;i++)
      os<<"\t";
    os<<(*itreenode)->getVertex()<<endl;
    itreenode++;
  }
  return os;
}



void TreeGraph::delSubTree(int vertex) 
{
  assert (vertex != 0);
    
  vector<TreeNode*>::iterator itreenode;
  vector<NodeList*>::iterator i_node_table;
  itreenode = _treenodes.begin();
  i_node_table = _outNodeTable.begin();
  int i;
  int nb_desc = getNbDesc(vertex);
  _depth=0;
  _degree=0;
  for (i=0;i<father(vertex);i++){
    if ((*itreenode)->depth()>_depth)
      _depth = (*itreenode)->depth();
    if ((*i_node_table)->size()>_degree)
      _degree = (*i_node_table)->size();
    itreenode++;
    i_node_table++;
  }
  /* on enlève vertex de la liste des fils de son père */
  removeSon((*itreenode)->getNumber(),vertex,nb_desc+1);
  while (i<vertex){
    itreenode++;
    i_node_table++;
    i++;
  }
  /* on efface vertex et tous ses descendants */
  //  cout<<"Nbre de descendants a effacer : "<<nb_desc<<endl;
  vector<TreeNode*>::iterator etreenode = itreenode;
  vector<NodeList*>::iterator e_node_table = i_node_table;
  for (i = 0; i<=nb_desc;i++){
    //    cout<<i<<endl;
    etreenode++;
    e_node_table++;
  }
  /* on renumérote le graphe */
  vector<TreeNode*>::iterator ntreenode = etreenode;
  vector<NodeList*>::iterator n_node_table = e_node_table;
  while (ntreenode != _treenodes.end()){
    (*ntreenode)->putNumber((*ntreenode)->getNumber()-nb_desc-1);
    vector<int>::iterator i_node_list;
    i_node_list = (*n_node_table)->begin();
    while (i_node_list != (*n_node_table)->end()){
      (*i_node_list) = (*i_node_list)-nb_desc-1;
      i_node_list++;
    }
    if ((*ntreenode)->depth()>_depth)
      _depth = (*ntreenode)->depth();
    if ((*i_node_table)->size()>_degree)
      _degree = (*n_node_table)->size();
    ntreenode++;
    n_node_table++;
  }
  _treenodes.erase(itreenode,etreenode);
  _outNodeTable.erase(i_node_table,e_node_table);
  _nbVertex -= nb_desc+1;
//   cout<<"Taille de _treenodes : "<<_treenodes.size()<<endl;
//   print();
}

void TreeGraph::addSubTree(TreeGraph* T,int root, int insertion_point) 
{
  /* on récupère dans les deux listes les points d'insertion */	     
  vector<TreeNode*>::iterator begin;
  vector<NodeList*>::iterator begin_node_table;
  begin = _treenodes.begin();
  begin_node_table = _outNodeTable.begin();
  for (int i=0;i<insertion_point;i++){
    begin++;
    begin_node_table++;
  }
  int new_number = insertion_point+1;
  int added_nb_vertex = T->getNbDesc(root)+1;
  int decalage = added_nb_vertex ;
  int insertion_depth = (getTreeNode(insertion_point))->depth();
  int root_depth = (T->getTreeNode(root))->depth();

  shiftSonNumber(father(insertion_point),insertion_point,decalage);
  // cout<<"Numero du noeud insere = "<<new_number<<endl;

  /*on rajoute un fils à insertion_point*/
  (*begin_node_table)->push_back(new_number);
  /* et on décale l'index de tous les fils d'insertion point */
  if (insertion_point !=0) /*ie on insère pas à la racine */
    shiftSonNumber(insertion_point,child(insertion_point,1),decalage);

  begin++; /* l'insertion se fait juste avant l'itérateur */
  begin_node_table++;
  int insert_decalage = new_number - root;
  for (int j=0;j<added_nb_vertex;j++){
    
    TreeNode* tnode = new TreeNode(*(T->getTreeNode(root+j)));
    tnode->putNumber(new_number);
    tnode->putDepth(tnode->depth()-root_depth+insertion_depth+1);
    tnode->putFather(tnode->father()+insert_decalage);
    _treenodes.insert(begin,tnode);

    /* On recupère la liste des fils et on leur met le bon index */
    NodeList* new_sonlist = new NodeList();
    NodeList sonlist =*(T->sons(root+j));
    vector<int>::iterator it_son_list;
    it_son_list = sonlist.begin();
    while (it_son_list != sonlist.end()){
      new_sonlist->push_back((*it_son_list) + insert_decalage);
      it_son_list++;
    }
    if (_degree<new_sonlist->size())
      _degree = new_sonlist->size();
    _outNodeTable.insert(begin_node_table,new_sonlist);
    new_number++;

/* On renumérote les derniers noeuds */ 
  }
  while (begin !=  _treenodes.end()){
    (*begin)->putNumber((*begin)->getNumber()+decalage);
    (*begin)->putFather((*begin)->father()+decalage);
    vector<int>::iterator i_node_list;
    i_node_list = (*begin_node_table)->begin();
    while (i_node_list != (*begin_node_table)->end()){
      (*i_node_list) = (*i_node_list)+decalage;
      i_node_list++;
    }
    begin++;
    begin_node_table++;
  }
  _nbVertex += added_nb_vertex;

}

int TreeGraph::minimalClass(int node) const{
  int mc = 10000;
  int nbSon = getNbChild(node);
  if(nbSon != 0){
    for (int i=1;i<= nbSon;i++){
      int son = child(node,i);
      if(_class[son] < mc)
	mc = _class[son];
    }
  }
  else
    mc = -1;
  return mc;
}


int TreeGraph::rightClasse(int clas){
  int rc,found,nbbrothers,i;
  vector<NodeList*>::const_iterator begin;  
  rc  = clas;
  found = 0;
  
  begin = _classNodeTable.begin();
  for (i=0;i<clas;i++)
    begin++;
  
  int node = *((*begin)->begin());
  
  nbbrothers = getNbChild(father(node));
  for (i=1;i<= nbbrothers;i++){
    int son = child(father(node),i);
    if(_class[son]> rc && found == 0){
      rc = _class[son];
      found =1;
    }
  }
  if (rc == clas)
    rc = -1;
  return rc;
}

NodeList* TreeGraph::getRootsInForestClass(int clas){
  NodeList* nlist = new NodeList();
  vector<NodeList*>::const_iterator begin;
  vector<int>::const_iterator begin_node;
  int current = clas;
  while(current != -1){
    begin = _classNodeTable.begin();
    for (int i=0;i<current;i++){
      begin++;
    }    
    
    begin_node = (*begin)->begin();
    while(begin_node !=  (*begin)->end()){
      nlist->push_back(*begin_node);
      begin_node++;
    }
    current = rightClasse(current);
  }  
  return nlist;
}

NodeList* TreeGraph::getRootsInClass(int clas){

  NodeList * nlist = new NodeList();
  vector<NodeList*>::const_iterator begin;
  vector<int>::const_iterator begin_node;
  
  begin = _classNodeTable.begin();
  for (int i=0;i<clas;i++){
    begin++;
  }    
  
  begin_node = (*begin)->begin();
  while(begin_node !=  (*begin)->end()){
    nlist->push_back(*begin_node);
    begin_node++;
  }
  
  return nlist;
}

int TreeGraph::getNbRootsInClass(int clas){
  int nb = 0;
  if(clas ==-1)
    nb = 0;
  else{
    nb = _classNodeTable[clas]->size();
  }
  return nb;
}

