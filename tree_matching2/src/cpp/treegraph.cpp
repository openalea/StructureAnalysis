/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       TreeMatching : Comparison of Tree Structures
 *
 *       Copyright 1995-2009 UMR LaBRI
 *
 *       File author(s): P.ferraro (pascal.ferraro@labri.fr)
 *
 *       $Source: /usr/cvsmaster/AMAPmod/src/TreeMatching/treegraph.cpp,v $
 *       $Id: treegraph.cpp 3258 2007-06-06 13:18:26Z dufourko $
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



TreeGraph::TreeGraph( )
{
  _nbVertex = 0;
  _rootId = 0;
  _degree = 0;
  _depth = 0;
  //_treenodes = 0;
}


//Destructor
TreeGraph::~TreeGraph()
{
  _treenodes.clear();
}

void TreeGraph::addNode(int id, int father)
{
  TreeNodePtr tree_node(new TreeNode(id,father));
  if (father == -1){
    _rootId = tree_node->getId();
    tree_node->putDepth(1);
    _depth = 1;
  }
  else{
    _treenodes[father]->addChild(id);
    int nb_child = _treenodes[father]->getChildNumber();
    if ( nb_child >_degree)
      _degree++;
    int depth = (_treenodes[father])->depth()+1;
    tree_node->putDepth(depth);
    if (depth > _depth)
      _depth = depth ;
  }
  _treenodes.push_back(tree_node);
}


void TreeGraph::addNode(TreeNodePtr tree_node)
{
  int father = tree_node->father();
  if (father == -1){
    _rootId = tree_node->getId();
    tree_node->putDepth(1);
    _depth = 1;
  }    
  else{
    _treenodes[father]->addChild(tree_node->getId());
    if (_treenodes[father]->getChildNumber() >_degree)
      _degree++;
    int depth = _treenodes[father]->depth()+1;
    tree_node->putDepth(depth);
    if (depth > _depth)
      _depth = depth ;
  }
  _treenodes.push_back(tree_node);
}

TreeNodePtr TreeGraph::getTreeNode(int vertex) 
{
  assert((vertex>=0)&&(vertex<_treenodes.size()));
  return _treenodes[vertex];
}

void TreeGraph::addValue(int index, DistanceType new_value){
  assert((index>=0)&&(index<_treenodes.size()));
  _treenodes[index]->addValue(new_value);
}

DistanceType TreeGraph::getValue(int index, int index_value) const {
  assert((index>=0)&&(index<_treenodes.size()));
  assert((index_value>=0)&&(index_value<_treenodes[index]->getValueSize()));
  return _treenodes[index]->getValue(index_value);
}


int TreeGraph::isLeaf(int node) const
{
  return((getNbChild(node)==0));
}


int TreeGraph::father(int node) const
{
  return (_treenodes[node])->father();
}

int TreeGraph::child(const int node,const int child_number) const
{
  if (child_number > getNbChild(node))
    return(-1);
  else{
    //   _treenodes[node].print();
    return (_treenodes[node]->getChild(child_number));
  }
}


const NodeList TreeGraph::childList(const int vertex) const
{
  return (_treenodes[vertex]->getChildList());
}

int TreeGraph::getNbChild(int vertex) const
{
  if (vertex==-1)
    return(0);

  else
    return(childList(vertex).size());
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

int TreeGraph::getNbVertex() const
{
  return(_treenodes.size());
}


TreeNodePtr TreeGraph::getNode(int i)
{
  return getTreeNode(i);
}


int TreeGraph::isNull()
{
  return (getNbVertex() == 0);
}

void TreeGraph::print() const
{
  cerr<<"Degre = "<<_degree<<endl;
  cerr<<"Vertex Number = "<<getNbVertex()<<endl;
  for (int i = 0; i < getNbVertex() ; i++)
    _treenodes[i]->print();
}

bool TreeGraph::mtg_write(  char *path ) 
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
  // print hearder of mtg file
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
  os<<"VId \t ALPHA"<<endl;
  os<<endl;
  os<<"MTG:"<<endl;
  os<<"ENTITY-CODE";
  int i;
  for (i=0;i<_depth;i++)
    os<<"\t";
  os<<"VId"<<endl;

  // print nodes
  TreeNodeList::const_iterator itreenode = _treenodes.begin();
  os<<"P1/";
  os<<"E"<< ((*itreenode)->getId());
  for (i=0;i<_depth-(*itreenode)->depth()+1;i++)
    os<<"\t";
  os<<"XX"<<endl;
  itreenode++;
  while (itreenode != _treenodes.end()){
    for (i=0;i<(*itreenode)->depth()-1;i++)
      os<<"\t";
    os<<"+E"<<(*itreenode)->getId();
    for (i=0;i<_depth-(*itreenode)->depth()+1;i++)
      os<<"\t";
    os<<"XX"<<endl;
    itreenode++;
  }
  return os;
}

