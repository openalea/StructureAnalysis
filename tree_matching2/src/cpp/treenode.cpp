/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       TreeMatching : Comparison of Tree Structures
 *
 *       Copyright 1995-2009 UMR LaBRI
 *
 *       File author(s): P.ferraro (pascal.ferraro@labri.fr)
 * 
 *       $Source$
 *       $Id: treenode.cpp 3258 2007-06-06 13:18:26Z dufourko $
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


#include "treenode.h"
using namespace std;


TreeNode::TreeNode(int id, int father)
{
  _id   = id;
  _father   = father;
  _values   = ValueVector();
  _depth = 0;
}

TreeNode::~TreeNode()
{

}


void TreeNode::print() const
{
  cout<<"ID \t : "<<_id<<endl;
  cout<<"FATHER \t : "<<_father<<endl;
  cout<<"DEPTH \t : "<<_depth<<endl;
  cout<<"CHILD LIST \t : [ ";
  for (int i=0;i<getChildNumber()-1;i++)
    cout<<getChild(i)<<" , ";
  if (getChildNumber()>0)
    cout<<getChild(getChildNumber()-1)<<" ] "<<endl;
  else
    cout<<"]"<<endl;

   for (int i=0;i<_values.size();i++)
     cout<<"ARG["<<i<<"] \t : "<<getValue(i)<<endl;
}


int TreeNode::getValueSize() const
{
  return(_values.size());
}

DistanceType TreeNode::getValue(int index) const
{
  if (_values.size())
  {
    assert((index>=0)&&(index<_values.size()));
    return(_values[index]);
  }
  else
  {
    return(DIST_UNDEF);
  }
}

void TreeNode::addValue(DistanceType new_value)
{
  _values.push_back(new_value);
}

void TreeNode::putValue(int index, DistanceType new_value)
{
  assert((index>=0)&&(index<_values.size()));
  if (_values.size()) 
    _values[index]=new_value;
  
}

void TreeNode::addChild(int child_id){
  _childList.push_back(child_id);
}

void TreeNode::setChildList(vector<int> child_list){
  _childList = child_list;
}

vector<int> TreeNode::getChildList() const{
  return _childList;
}

int TreeNode::getChild(int id) const{
  assert ((id>=0) && (id<_childList.size()));
  return _childList[id];
}

int TreeNode::getChildNumber() const{
  return _childList.size();
}

